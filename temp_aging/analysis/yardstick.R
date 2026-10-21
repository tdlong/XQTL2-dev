# yardstick.R — marker-based heritability of each trait: how much of the
# phenotypic variance founder haplotypes account for, genome-wide.
#
#   bash temp_aging/analysis/run_numbers.sh yardstick
#
# DATA: the same 10-replicate scan as every other number here
# (AGE_SY_4scan_no89, replicates 8 and 9 dropped), on the same 268 1 cM tiles
# as scan_resolution.R and h2_threshold.R. Nothing about the hold-out analysis
# enters the estimate; loo_u.txt.gz is read only for its per-replicate values,
# which give a jackknife standard error over replicates.
#
# THE ESTIMATE. Summing the unbiased tile h2 over the genome counts one locus at
# every tile it is linked to. The number of tiles a locus in tile t is counted
# in is its LD score ell_t, measured from drift in ld_local.R. Marker-based
# heritability is sum_t h2_t / ell_t -- LD score regression's correction with
# the LD score measured here -- and is the counterpart of SNP-based heritability
# at the resolution of founder haplotype blocks.

suppressMessages(library(tidyverse))

SCAN <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
U    <- "temp_aging/numbers/loo_u.txt.gz"
TS1  <- "temp_aging/analysis/TableS1.tsv"
TSV  <- "temp_aging/numbers/yardstick.tsv"     # the table make_LOO_fig.R reads
PEAK <- c(chr = "chr3L", cM = 28.2)             # Table S1
TRTS <- c("SY10_F", "SY10_M", "SY20_F", "SY20_M")

d <- read.table(SCAN, header = TRUE, sep = "\t") %>% as_tibble() %>%
  mutate(trt = paste0(sugar, "_", sex))
HET <- read.table("pipeline/helpfiles/het_bounds.txt", header = TRUE, comment.char = "#")

# ── tiles, exactly as h2_threshold.R / scan_resolution.R ─────────────────────
w <- d %>% left_join(HET, by = "chr") %>% group_by(chr) %>%
  mutate(tile = floor((cM - min(cM)) / 1),
         is_eu = pos/1e6 >= eu_start & pos/1e6 <= eu_end) %>% ungroup()
keep <- w %>% group_by(chr, tile) %>% summarise(f = mean(is_eu), .groups = "drop") %>%
  filter(f > 0) %>% select(chr, tile)
w <- w %>% semi_join(keep, by = c("chr", "tile")) %>% filter(is_eu)
tiles <- w %>% group_by(trt, chr, tile) %>%
  summarise(h2 = mean(H2), wald = max(Wald_log10p), .groups = "drop")
cat("tiles:", n_distinct(paste(tiles$chr, tiles$tile)), "of 1 cM\n")

# ── the LD score of each tile, from drift (ld_local.R) ───────────────────────
# ell_t = number of 1 cM tiles a locus in tile t is counted in, measured from
# the cage's movement between draws, not from any selection signal. The 3L 9.30
# footprint (3.9-4.2) was used here before; drift shows that locus does NOT sit
# in a long-block region (ell = 1.1 there), so its width is not LD.
LD <- read.table("temp_aging/numbers/ld_local.tsv", header = TRUE, sep = "\t") %>% as_tibble()
tiles <- tiles %>% left_join(w %>% group_by(chr, tile) %>% summarise(cMmid = mean(cM), .groups = "drop"),
                             by = c("chr", "tile")) %>%
  mutate(ltile = floor(cMmid)) %>% left_join(LD %>% select(chr, ltile = tile, ell), by = c("chr", "ltile")) %>%
  group_by(chr) %>% mutate(ell = ifelse(is.na(ell), median(ell, na.rm = TRUE), ell)) %>% ungroup()
cat(sprintf("\nLD score per tile: median %.2f, 10-90%% %.2f-%.2f (ld_local.R)\n",
            median(tiles$ell), quantile(tiles$ell, .1), quantile(tiles$ell, .9)))

# ── the yardstick, and where it comes from ───────────────────────────────────
ys <- tiles %>% group_by(trt) %>% summarise(sum_tiles = sum(h2), h2 = sum(h2 / ell), .groups = "drop")
TPL <- 1   # kept for the band table below, which divides each tile by its own ell

cat("\nmarker-based heritability (%): sum over the genome of tile h2 / ell_t\n\n")
band <- tiles %>% mutate(band = cut(wald, c(-Inf, 2, 7.5, 15, Inf),
                                    c("LWP < 2", "2 to 7.5", "7.5 to 15", "> 15"))) %>%
  group_by(trt, band) %>% summarise(n = n(), h2 = sum(h2 / ell), .groups = "drop")
tab <- band %>% pivot_wider(names_from = trt, values_from = c(n, h2)) %>%
  select(band, starts_with("n_"), starts_with("h2_"))
print(as.data.frame(bind_rows(
  tab %>% mutate(across(where(is.double), ~round(.x, 2))),
  tibble(band = "whole genome") %>%
    bind_cols(as_tibble(as.list(setNames(tiles %>% count(trt) %>% pull(n), paste0("n_", TRTS))))) %>%
    bind_cols(as_tibble(as.list(setNames(round(ys$h2, 2), paste0("h2_", TRTS))))))),
  row.names = FALSE)
cat("\n  n = tiles in the band (by the tile's largest Wald in that treatment);\n",
    " h2 = that band's share of the marker-based heritability.\n")

# ── jackknife SE over replicates, from the per-replicate values in loo_u ─────
u <- read.table(U, header = TRUE, sep = "\t") %>% as_tibble() %>%
  transmute(trt = paste0(sugar, "_", sex), chr, tile, founder, k, rep = holdout, u = u_te)
u <- u %>% left_join(LD %>% select(chr, tile, ell), by = c("chr", "tile")) %>%
  group_by(chr) %>% mutate(ell = ifelse(is.na(ell), median(ell, na.rm = TRUE), ell)) %>% ungroup()
full <- function(z) z %>% group_by(trt, chr, tile, founder, k, ell) %>%
  summarise(m = mean(u), v = var(u), n = n(), .groups = "drop") %>%
  mutate(h = 100 * k * (m^2 - v/n) / ell) %>% group_by(trt) %>%
  summarise(h2 = sum(h), .groups = "drop")
jk <- u %>% group_split(trt) %>% map_dfr(function(z) {
  th <- map_dbl(unique(z$rep), ~ full(z %>% filter(rep != .x))$h2)
  tibble(trt = z$trt[1], from_u = full(z)$h2,
         se = sqrt((length(th) - 1) / length(th) * sum((th - mean(th))^2))) })
ys <- ys %>% left_join(jk, by = "trt")
cat("\nstandard error: jackknife over the ten replicates\n\n")
print(as.data.frame(ys %>% transmute(trt, h2 = round(h2, 1), se = round(se, 1),
  from_loo_u = round(from_u, 1))), row.names = FALSE)
cat("\n  from_loo_u recomputes the estimator from per-replicate values on\n",
    " loo_genome.R's 273 tiles (heterochromatin not dropped); it is the\n",
    " cross-check that the two files agree, and the source of the SE.\n")

# ── against Table S1 ─────────────────────────────────────────────────────────
s1 <- read.table(TS1, header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE) %>%
  as_tibble() %>% filter(!is.na(`Position (Mb)`)) %>%
  mutate(set = if_else(row_number() <= 8, "LWP > 15 (8 peaks)", "7.5 to 15 (10 peaks)"))
s1l <- s1 %>% select(set, starts_with("h2 ")) %>%
  pivot_longer(-set, names_to = "trt", values_to = "h2") %>%
  mutate(trt = sub("h2 (SY\\d+) ([FM])", "\\1_\\2", trt))
cat("\nthe peaks of Table S1: sum of their h2, and as a fraction of the yardstick\n\n")
print(as.data.frame(bind_rows(s1l %>% group_by(set, trt) %>% summarise(h2 = sum(h2), .groups = "drop"),
  s1l %>% group_by(trt) %>% summarise(h2 = sum(h2), .groups = "drop") %>% mutate(set = "all 18 peaks")) %>%
  left_join(ys %>% select(trt, total = h2), by = "trt") %>%
  mutate(frac = round(h2 / total, 2), h2 = round(h2, 2)) %>% select(-total) %>%
  pivot_wider(names_from = trt, values_from = c(h2, frac))), row.names = FALSE)
cat("\n  Peak h2 is at the peak window and chosen in-sample, so these are upper\n",
    " bounds on what the peaks explain; the held-out counterpart is Fig. Yb.\n")

write.table(ys %>% transmute(trt, h2, se, sum_tiles), TSV,
            row.names = FALSE, quote = FALSE, sep = "\t")
cat("\nWrote:", TSV, "\n")
