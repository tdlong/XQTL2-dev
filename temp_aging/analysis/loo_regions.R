# loo_regions.R — fraction of marker-based heritability a predictor explains in
# the held-out replicate, counted so that no LD correction is needed.
#
#   bash temp_aging/analysis/run_numbers.sh loo_regions
#
# THE PROBLEM IT AVOIDS. Summing tile h2 counts a locus once per tile its
# signal reaches, and how many that is (the LD score) is contested: drift says
# ~1.3 tiles, the signal's own footprint at a big peak ~4. Any comparison of a
# predictor on ISOLATED tiles against the genome-wide SUM inherits that factor.
#
# THE FIX. Let a predictor's marker be a REGION: the chosen tile plus every tile
# within +-W cM of it. Its explained h2 and the genome-wide sum then count tiles
# the same way and the over-count cancels in the ratio. At W = all tiles the
# predictor is an all-marker PRS and the ratio is exact.
#
# PER FOLD (loo_genome.R's loo_tiles.txt.gz): centres are the best tiles on the
# nine, 5 cM apart, as before. For a region set S,
#   cross = sum_S 100 k sum_f e_bar_f e_f(r) / p_f      (h2 in the held-out replicate)
#   nine  = sum_S 100 k sum_f e_bar_f^2 / p_f           (the predictor's variance)
#   explained = mean(cross)^2 / mean(nine) over the ten folds  (cov^2 / var)
# divided by the genome-wide sum of the full ten-replicate tile h2.

suppressMessages(library(tidyverse))
TIL <- "temp_aging/numbers/loo_tiles.txt.gz"
UF  <- "temp_aging/numbers/loo_u.txt.gz"
TSV <- "temp_aging/numbers/loo_regions.tsv"
WS  <- c(0, 1, 2, 3, 5); NS <- 1:30; CLUMP <- 5

T <- read.table(TIL, header = TRUE, sep = "\t") %>% as_tibble() %>% mutate(trt = paste0(sugar, "_", sex))
U <- read.table(UF, header = TRUE, sep = "\t") %>% as_tibble() %>% mutate(trt = paste0(sugar, "_", sex))
tot <- U %>% group_by(trt, chr, tile, founder, k) %>%
  summarise(m = mean(u_te), v = var(u_te), n = n(), .groups = "drop") %>%
  mutate(h = 100 * k * (m^2 - v / n)) %>% group_by(trt) %>% summarise(sum_full = sum(h), .groups = "drop")
cat("tiles:", n_distinct(paste(T$chr, T$tile)), " folds:", n_distinct(T$holdout), "\n")
cat("\ngenome-wide sum of tile h2, full ten replicates (the denominator; % x LD score)\n")
print(as.data.frame(tot %>% mutate(sum_full = round(sum_full, 1))), row.names = FALSE)

clump <- function(z, win, n) { out <- integer(0); z <- z[order(-z$nine), ]
  while (nrow(z) && length(out) < n) { out <- c(out, z$id[1])
    z <- z[!(z$chr == z$chr[1] & abs(z$cM - z$cM[1]) < win), ] }; out }
folds <- T %>% group_split(trt, holdout) %>% map_dfr(function(z) {
  z <- z %>% mutate(id = row_number()); cen <- clump(z, CLUMP, max(NS))
  dmin <- function(c) map_dbl(seq_len(nrow(z)), function(i) {        # distance of each tile to nearest centre
    d <- abs(z$cM[i] - c$cM); d[z$chr[i] != c$chr] <- Inf; min(d) })
  per_N <- map(NS, function(N) dmin(z[cen[seq_len(min(N, length(cen)))], ]))
  bind_rows(
    map_dfr(WS, function(W) map_dfr(seq_along(NS), function(j) { inS <- per_N[[j]] <= W
      tibble(W = W, N = NS[j], cross = sum(z$cross[inS]), nine = sum(z$nine[inS]), tiles = sum(inS)) })),
    tibble(W = Inf, N = Inf, cross = sum(z$cross), nine = sum(z$nine), tiles = nrow(z))) %>%
    mutate(trt = z$trt[1], holdout = z$holdout[1]) })
# jackknife over folds of the ratio-of-means estimator itself
jk_ratio <- function(cross, nine) { n <- length(cross)
  th <- map_dbl(seq_len(n), ~ mean(cross[-.x])^2 / mean(nine[-.x]))
  sqrt((n - 1) / n * sum((th - mean(th))^2)) }
out <- folds %>% left_join(tot, by = "trt") %>% group_by(trt, W, N) %>%
  summarise(explained = mean(cross)^2 / mean(nine), tiles = mean(tiles),
            frac = explained / sum_full[1], se = jk_ratio(cross, nine) / sum_full[1], .groups = "drop")

cat("\nfraction of marker-based heritability explained in the held-out replicate\nby the top N regions of +-W cM (W = 0: single tiles; Inf: every tile)\n\n")
print(as.data.frame(out %>% filter(N %in% c(1, 5, 10, 20, 30, Inf)) %>% mutate(frac = round(frac, 2)) %>%
  select(trt, W, N, frac) %>% pivot_wider(names_from = trt, values_from = frac) %>% arrange(W, N)), row.names = FALSE)
cat("\ntiles in the predictor, mean over folds (SY10_F)\n\n")
print(as.data.frame(out %>% filter(trt == "SY10_F", N %in% c(1, 5, 10, 20, 30, Inf)) %>%
  transmute(W, N, tiles = round(tiles)) %>% pivot_wider(names_from = N, values_from = tiles)), row.names = FALSE)
cat("\n  The W = Inf row is the ceiling: an all-marker predictor built on nine\n",
    " replicates explains this fraction of what the tenth carries. At W = 0 the\n",
    " thirty markers are 30 of 273 tiles; at W = 3 they cover ~160.\n")
write.table(out, TSV, row.names = FALSE, quote = FALSE, sep = "\t")
cat("\nWrote:", TSV, "\n")
