# ld_drift.R — how far a founder haplotype block extends, measured from drift.
#
#   bash temp_aging/analysis/run_numbers.sh ld_drift
#
# WHY. Summing tile h2 over the genome counts a locus once per tile it is
# linked to. The number of tiles is the LD score, sum over distance of g(d)^2,
# with g(d) the correlation between a locus and a tile d cM away. It can be read
# off a big peak's footprint, but that assumes one footprint for the genome.
# Here it is measured independently of any selection signal. There is ONE cage,
# and the ten replicates are cohorts drawn from it at different times (two may
# share a date). Its founder frequencies move between draws -- drift, and
# whatever selection the cage itself is under -- and a move at one window is
# carried to linked windows with the same g(d). A replicate's deviation from
# the mean over the ten draws is the cage's state at that time. It is measured
# twice, by the replicate's control pool and its selected pool, drawn from the
# same cohort with independent sampling noise, so the cross product of their
# deviations is noise-free at every lag, lag 0 included. C(d) is then the
# autocorrelation of g, deconvolved to g.
#
# DATA. AGE_SY_zoom_means: per-replicate, per-pool founder frequencies in 5 kb
# windows over seven regions of 1.2 Mb around peaks (replicates 8, 9 dropped).
# Controls are one pool per replicate per sex, listed under both sugars.

suppressMessages(library(tidyverse))

Z <- read.table("process/AGE_SY/AGE_SY_zoom_means.txt.gz", header = TRUE, sep = "\t") %>%
  as_tibble() %>% filter(!REP %in% c(8, 9))
cm <- read.table("process/AGE_SY/AGE_SY_4scan_no89.txt.gz", header = TRUE, sep = "\t") %>%
  distinct(chr, pos, cM)
Z <- Z %>% inner_join(cm, by = c("chr", "pos"))
cat("regions:", n_distinct(Z$locus), " replicates:", n_distinct(Z$REP), " windows per region:",
    n_distinct(Z$pos) / n_distinct(Z$locus), "\n")

C <- Z %>% filter(TRT == "C", sugar == "SY10") %>% group_by(locus, sex, founder, pos) %>%
  mutate(dc = freq - mean(freq)) %>% ungroup() %>% select(locus, sex, founder, REP, pos, cM, dc)
S <- Z %>% filter(TRT == "Z") %>% group_by(locus, sugar, sex, founder, pos) %>%
  mutate(dz = freq - mean(freq)) %>% ungroup() %>% select(locus, sugar, sex, founder, REP, pos, dz)
D <- C %>% inner_join(S, by = c("locus", "sex", "founder", "REP", "pos")) %>%
  arrange(locus, sugar, sex, founder, REP, cM)

xlag <- function(D, L) {
  w <- D %>% group_by(locus, sugar, sex, founder, REP) %>%
    mutate(j = findInterval(cM + L, cM),
           ok = j >= 1 & j <= n() & abs(cM[pmax(j, 1)] - (cM + L)) < 0.05,
           dz2 = dz[pmax(j, 1)], dc2 = dc[pmax(j, 1)]) %>% ungroup() %>% filter(ok)
  tibble(lag = L, s = mean(c(w$dc * w$dz2, w$dz * w$dc2)), n = nrow(w))
}
LAGS <- seq(0, 4, by = 0.5)
A <- map_dfr(LAGS, ~xlag(D, .x)) %>% mutate(rho = s / s[lag == 0])
cat("\nC(d): autocorrelation along the genome of the cage's state at each draw, all regions pooled\n\n")
print(as.data.frame(A %>% mutate(rho = round(rho, 3), s = signif(s, 3))), row.names = FALSE)

cat("\nby region (lag in cM)\n\n")
B <- D %>% group_split(locus) %>% map_dfr(function(z)
  map_dfr(0:4, ~xlag(z, .x)) %>% mutate(locus = z$locus[1], rho = round(s / s[lag == 0], 2)) %>%
    select(locus, lag, rho)) %>%
  pivot_wider(names_from = lag, values_from = rho, names_prefix = "d")
print(as.data.frame(B), row.names = FALSE)
cat("\n  NaN: the region is too short for that lag. 3L 9.31 is wider than the rest.\n")

# deconvolve: C = (g*g)/sum g^2, so the spectrum of C is |ghat|^2
deconv <- function(rho, step) {
  N <- 512; L <- length(rho); c0 <- numeric(N)
  c0[1:L] <- rho; c0[N:(N - L + 2)] <- rho[2:L]
  sp <- Re(fft(c0)); sp[sp < 0] <- 0
  g <- Re(fft(sqrt(sp), inverse = TRUE)) / N; g <- g / g[1]
  sum(g^2) * step                                  # LD score in 1 cM tiles
}
ld_all <- deconv(pmax(A$rho, 0), 0.5)
ld_3L  <- deconv(pmax(B %>% filter(locus == "chr3L:9.31") %>% select(starts_with("d")) %>% unlist(), 0), 1)
ld_not <- deconv(pmax(map_dbl(0:4, ~xlag(D %>% filter(locus != "chr3L:9.31"), .x)$s) /
                      xlag(D %>% filter(locus != "chr3L:9.31"), 0)$s, 0), 1)
cat(sprintf("\nLD score (tiles a locus is counted in) from the cage's movement between draws: all regions %.1f, 3L 9.31 %.1f, the other six %.1f\n",
            ld_all, ld_3L, ld_not))
cat("For comparison: the 3L 9.30 footprint of the selection signal gives 3.9-4.2 (yardstick.R);\n",
    "the cross-fold autocorrelation of the signal genome-wide, deconvolved the same way, 2.2-2.3.\n")
