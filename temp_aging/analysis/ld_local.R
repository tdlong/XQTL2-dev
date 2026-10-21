# ld_local.R — the LD score of every 1 cM tile, from the cage's movement
# between draws, genome-wide. Extends ld_drift.R (seven regions) to the genome.
#
#   bash temp_aging/analysis/run_numbers.sh ld_local
#
# WHAT. ell_t = integral over distance of g(d)^2, with g(d) the correlation
# between a locus in tile t and a point d cM away: the number of 1 cM tiles a
# locus in t is counted in when tile h2 is summed. Marker-based heritability is
# sum_t h2_t / ell_t (yardstick.R).
#
# HOW. One cage, ten draws. A replicate's deviation from the ten-draw mean at a
# window is the cage's state at that draw; it is carried to linked windows with
# g(d). The control and selected pools of a replicate sample the same cohort
# with independent noise, so dc(x) * dz(x + d) averaged over founders,
# replicates, sexes and sugars is the drift autocovariance at lag d, noise-free
# at every lag including 0. Normalised, that is C(d), the autocorrelation of g;
# g is recovered by deconvolution (spectrum of C = |spectrum of g|^2) and
# ell = sum g^2 * step. Done once genome-wide, and once per tile using the
# windows within +-HALF cM of it.
#
# DATA. drift_tiles.txt.gz: per-replicate pool founder frequencies on a 0.25 cM
# grid (drift_tiles.R, HPC3).

suppressMessages(library(tidyverse))

IN   <- "temp_aging/numbers/drift_tiles.txt.gz"
TSV  <- "temp_aging/numbers/ld_local.tsv"
STEP <- 0.25; MAXLAG <- 10; HALF <- 5

D <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble()
cat("grid points:", n_distinct(paste(D$chr, D$cM)), " replicates:", n_distinct(D$REP), "\n")

# deviations from the ten-draw mean; controls are shared by the two sugars
dc <- D %>% filter(pool == "C") %>% group_by(chr, cM, sex, founder) %>%
  mutate(dev = freq - mean(freq)) %>% ungroup() %>% select(chr, cM, sex, founder, REP, dc = dev)
dz <- D %>% filter(pool == "Z") %>% group_by(chr, cM, sex, sugar, founder) %>%
  mutate(dev = freq - mean(freq)) %>% ungroup() %>% select(chr, cM, sex, sugar, founder, REP, dz = dev)
P <- dz %>% inner_join(dc, by = c("chr", "cM", "sex", "founder", "REP"))   # one row per (point, sex, sugar, founder, rep)

# product sums per START point and lag, over everything else
key <- function(chr, cM) paste(chr, sprintf("%.2f", cM))
P <- P %>% mutate(k = key(chr, cM))
lags <- seq(0, MAXLAG, by = STEP)
prod_by_start <- map_dfr(lags, function(L) {
  Q <- P %>% mutate(k2 = key(chr, cM + L)) %>%
    inner_join(P %>% select(k2 = k, sex, sugar, founder, REP, dc2 = dc, dz2 = dz),
               by = c("k2", "sex", "sugar", "founder", "REP"))
  Q %>% group_by(chr, cM) %>% summarise(s = sum(dc * dz2 + dz * dc2) / 2, n = n(), .groups = "drop") %>%
    mutate(lag = L)
})

deconv <- function(rho, step) {                    # rho on a regular grid from lag 0
  z <- which(rho < 0); if (length(z)) rho[min(z):length(rho)] <- 0   # stop at the first zero crossing
  N <- 1024; L <- length(rho); c0 <- numeric(N)
  c0[1:L] <- rho; c0[N:(N - L + 2)] <- rho[2:L]
  sp <- Re(fft(c0)); sp[sp < 0] <- 0
  g <- Re(fft(sqrt(sp), inverse = TRUE)) / N; g <- g / g[1]
  sum(g^2) * step
}
# exponential block, g(d) = exp(-d/lambda): C(d) = (1 + d/lambda) exp(-d/lambda) and
# ell = integral g^2 = lambda. Fitted on lags 0.25-5 cM. (The deconvolution runs
# about 15% low on a synthetic exponential; the fit is exact on it.)
lam_fit <- function(rho, lags, lo = 0.25, hi = 5) {
  i <- lags >= lo & lags <= hi; r <- rho[i]; d <- lags[i]
  optimize(function(l) sum((r - (1 + d/l) * exp(-d/l))^2), c(0.2, 15))$minimum
}
rho_of <- function(S) { S <- S %>% group_by(lag) %>% summarise(s = sum(s), n = sum(n), .groups = "drop") %>%
  arrange(lag); S$s / S$s[1] }

cat("\ngenome-wide C(d), both sexes pooled\n\n")
rg <- rho_of(prod_by_start)
show <- c(0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10)
print(as.data.frame(tibble(lag_cM = show, C = round(rg[match(show, lags)], 3))), row.names = FALSE)
lg <- lam_fit(rg, lags)
cat(sprintf("\n  LD score genome-wide: deconvolved %.2f; exponential block lambda = %.2f cM (= ell) with\n  fitted C(1..4) = %s\n",
    deconv(rg, STEP), lg, paste(round((1 + (1:4)/lg) * exp(-(1:4)/lg), 3), collapse = " ")))

# per tile: products whose start lies within +-HALF cM of the tile centre
tiles <- prod_by_start %>% distinct(chr, cM) %>% mutate(tile = floor(cM)) %>% distinct(chr, tile) %>% arrange(chr, tile)
ell <- pmap_dfr(tiles, function(chr, tile) {
  S <- prod_by_start %>% filter(chr == !!chr, abs(cM - (tile + 0.5)) <= HALF)
  r <- rho_of(S)
  tibble(chr = chr, tile = tile, n_points = n_distinct(S$cM), ell_deconv = deconv(r, STEP), ell = lam_fit(r, lags))
})
cat("\nper-tile LD score ell_t = lambda_t, exponential fit (windows of +-", HALF, "cM)\n\n")
print(as.data.frame(ell %>% group_by(chr) %>% summarise(tiles = n(), median = round(median(ell), 2),
  q10 = round(quantile(ell, .1), 2), q90 = round(quantile(ell, .9), 2), max = round(max(ell), 2))), row.names = FALSE)
cat("\n  genome: median", round(median(ell$ell), 2), " range", paste(round(range(ell$ell), 2), collapse = "-"), "\n")

cat("\nell_t at the Table S1 peaks\n\n")
pk <- tribble(~peak, ~chr, ~cM, "3L 9.30", "chr3L", 28.2, "3R 8.63", "chr3R", 48.2, "3L 21.62", "chr3L", 46.9,
  "2L 11.95", "chr2L", 44.4, "2L 14.85", "chr2L", 50.0, "2R 13.88", "chr2R", 69.6, "2R 24.08", "chr2R", 103.3, "3L 7.43", "chr3L", 20.6)
print(as.data.frame(pk %>% mutate(tile = floor(cM)) %>% left_join(ell, by = c("chr", "tile")) %>%
  transmute(peak, ell = round(ell, 2), ell_deconv = round(ell_deconv, 2))), row.names = FALSE)

write.table(ell %>% select(chr, tile, ell, ell_deconv), TSV, row.names = FALSE, quote = FALSE, sep = "\t")
cat("\nWrote:", TSV, "\n")
