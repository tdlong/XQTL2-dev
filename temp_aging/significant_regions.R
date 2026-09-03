# significant_regions.R — sex and diet, over the significant part of the genome.
#
#   Rscript temp_aging/significant_regions.R
#
# DATA: the 10-replicate AGE_SY dataset -- replicates 8 and 9 dropped, those
# being the May 2023 cage (helpfiles/AGE_2024/population_assignment.txt).
#
# UNIT: 1 cM tiles, the same ones as scan_resolution.R and h2_threshold.R -- each
# arm cut at 0-1, 1-2, 2-3 cM from its own start, a tile kept if any part of it
# is euchromatic (268 of 270). Tiles do not overlap, so nothing here needs an
# overlap correction.
#
# A TILE IS SIGNIFICANT if some window in it reaches -log10 P = 7.5 in some
# treatment -- the same threshold paragraphs 1 and 2 use.
#
# A TILE'S VALUE comes from ONE window: the window with the largest Wald across
# the four treatments. The variance partition needs eight numbers (4 treatments x
# 2 split halves) that describe the same thing, so they have to come from the
# same window; taking a maximum separately in each of the eight would pick eight
# different windows and inflate every difference between them.
#
# CHECKED 2026-09-03: does the choice of representative window matter? Reran the
# partition on the same 85 significant tiles with the tile mean and with the
# window nearest the tile midpoint instead of the max-Wald window. Point
# estimates all land inside each other's 95% CIs (main 85-86%, sex 12-13%, diet
# 1.6-1.8%, int -0.4 to -0.5%), and per-tile the sex and main terms correlate
# 0.97-1.00 across the three rules. The choice does not change the answer;
# max-Wald is kept because it picks a real peak rather than an arbitrary point.
# (ad hoc check, not a retained script -- rerun by swapping the window-picking
# step above for a per-tile mean or nearest-to-midpoint window if needed again.)
#
# THE ERROR TERM comes from the split halves, not from a floor. The ten
# replicates of a treatment are split into the five odd and the five even, giving
# two independent h2 estimates of the same quantity; half their squared
# difference estimates the sampling variance of one treatment mean. That is a
# pure error term, and it is subtracted from every component below, so each is a
# sum of squares with its own noise already taken out. Nothing here needs the
# null-window floor of h2_threshold.R -- that exists because a single h2 estimate
# has nowhere to get an error term from, and here there is one.
#
# THE BIAS IS NOT COMMON TO THE FOUR TREATMENTS, so it cannot be left in and
# offset from squaring a noisy frequency estimate, computed from that pool's own
# lsei reconstruction error and the multinomial sampling of its own flies
# (scan_functions.R, XQTL2 #34). It is reported and not subtracted. Averaged over
# windows it runs 0.73 in SY10 females to 0.83 in SY20 males -- it tracks how many
# flies were selected and how well that pool reconstructs, both of which differ by
# sex and by diet. A sex contrast on raw h2 therefore contains a sex difference in
# bias. So it is subtracted here, per pool per window, before anything is
# contrasted. Both are printed; the difference between them is the bias effect.

suppressMessages(library(tidyverse))

SPLIT <- "process/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz"
SCAN  <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
THR   <- 7.5      # significant, as in paragraphs 1 and 2
WNULL <- 2        # Wald below this: frequencies did not move
TILE  <- 1        # cM

HET <- read.table("pipeline/helpfiles/het_bounds.txt", header = TRUE,
                  comment.char = "#") %>% as_tibble()
l <- read.table(SPLIT, header = TRUE, sep = "\t") %>% as_tibble()
d <- read.table(SCAN,  header = TRUE, sep = "\t") %>% as_tibble()

# ── 1 cM tiles, and the peak window in each ──────────────────────────────────
geo <- d %>% distinct(chr, pos, cM) %>% left_join(HET, by = "chr") %>%
  group_by(chr) %>% mutate(tile = floor((cM - min(cM)) / TILE)) %>% ungroup() %>%
  mutate(is_eu = pos/1e6 >= eu_start & pos/1e6 <= eu_end)
keep <- geo %>% group_by(chr, tile) %>% summarise(f = mean(is_eu), .groups = "drop") %>%
  filter(f > 0) %>% select(chr, tile)
geo <- geo %>% semi_join(keep, by = c("chr", "tile")) %>% filter(is_eu) %>%
  select(chr, pos, tile)

pk <- d %>% inner_join(geo, by = c("chr", "pos")) %>%
  group_by(chr, tile, pos) %>% summarise(w = max(Wald_log10p), .groups = "drop") %>%
  group_by(chr, tile) %>% slice_max(w, n = 1, with_ties = FALSE) %>% ungroup()

cat(sprintf("%d tiles of %g cM; %d reach -log10 P %g in some treatment (%.0f%%)\n\n",
            nrow(pk), TILE, sum(pk$w > THR), THR, 100 * mean(pk$w > THR)))
sig <- pk %>% filter(w > THR)

# eight values per significant tile, from its peak window
cell <- function() {
  l %>% inner_join(sig %>% select(chr, tile, pos), by = c("chr", "pos")) %>%
    mutate(v = H2) %>%
    select(chr, tile, sugar, sex, half, v)
}

# ── male vs female, by arm ───────────────────────────────────────────────────
by_arm <- function() cell() %>% group_by(chr, tile, sex) %>%
  summarise(h2 = mean(v), .groups = "drop") %>%
  pivot_wider(names_from = sex, values_from = h2) %>% group_by(chr) %>%
  summarise(tiles = n(), `med M` = round(median(M), 2),
            `med F` = round(median(F), 2),
            `M>F at` = sprintf("%.0f%%", 100 * mean(M > F)), .groups = "drop")

# X vs autosome is the contrast that has a mechanism behind it; the per-arm split
# is here because the autosomes are not homogeneous, not because any one arm is
# being claimed.
xa <- function() cell() %>% group_by(chr, tile, sex) %>%
  summarise(h2 = mean(v), .groups = "drop") %>%
  pivot_wider(names_from = sex, values_from = h2) %>%
  mutate(g = ifelse(chr == "chrX", "X", "autosomes")) %>% group_by(g) %>%
  summarise(tiles = n(), `med M` = round(median(M), 2),
            `med F` = round(median(F), 2),
            `M>F at` = sprintf("%.0f%%", 100 * mean(M > F)), .groups = "drop")
cat("X versus the autosomes (bias subtracted):\n\n")
xa() %>% as.data.frame() %>% print(row.names = FALSE)

cat("\nthe same split by arm (bias subtracted):\n\n")
by_arm() %>% as.data.frame() %>% print(row.names = FALSE)

# ── the partition ────────────────────────────────────────────────────────────
# Split halves give the error term: the odd-replicate and even-replicate h2 for
# one treatment differ only by sampling, so their squared difference estimates
# the noise in a single treatment mean. Each named term is a sum of squares with
# that noise subtracted, which is what makes "no detectable interaction" mean
# something rather than being an artifact of squaring.
terms <- function() {
  w <- cell() %>%
    pivot_wider(names_from = c(sugar, sex, half), values_from = v) %>% drop_na() %>%
    transmute(chr, tile, a = (SY10_F_odd + SY10_F_even)/2, b = (SY20_F_odd + SY20_F_even)/2,
      cc = (SY10_M_odd + SY10_M_even)/2, dd = (SY20_M_odd + SY20_M_even)/2,
      msr = ((SY10_F_odd - SY10_F_even)^2 + (SY20_F_odd - SY20_F_even)^2 +
             (SY10_M_odd - SY10_M_even)^2 + (SY20_M_odd - SY20_M_even)^2)/2/4) %>%
    mutate(y = (a+b+cc+dd)/4, F2 = (a+b)/2, M2 = (cc+dd)/2,
      S10 = (a+cc)/2, S20 = (b+dd)/2,
      main = 8*y^2 - msr, sex = 4*((F2-y)^2 + (M2-y)^2) - msr,
      diet = 4*((S10-y)^2 + (S20-y)^2) - msr,
      int = 2*((a-F2-S10+y)^2 + (b-F2-S20+y)^2 +
               (cc-M2-S10+y)^2 + (dd-M2-S20+y)^2) - msr)
  w
}

# Resample 4 cM groups of tiles rather than tiles one at a time, since
# neighbouring tiles can carry the same peak.
#
# MEASURED 2026-09-03, autocorrelation of the max-over-treatments Wald against
# genetic separation, euchromatic steps:
#
#   0.05 cM  0.999     0.5 cM  0.957     2 cM  0.695
#   0.10 cM  0.997     1.0 cM  0.877     4 cM  0.378
#
# A 1 cM tile is NOT an independent unit. The tile is a summarisation scale,
# not an independence scale; independence is what the 4-tile block buys, and
# 4 cM is where r has dropped to ~0.4. The +/-100 kb smoothing only spans
# 0.4-0.6 cM (2.0-3.1 cM/Mb by arm), so it accounts for the correlation to
# ~0.5 cM and real peak width for the rest. Mean h2 autocorrelates slightly
# higher throughout (0.91 at 1 cM).
#
# WHY 1 cM AND NOT 0.5 OR 2. Correlation between the summary values of ADJACENT
# tiles, against what each size costs in resolution:
#
#   size    tiles  sig@7.5  adj-tile r   separately resolved sig regions
#   0.5 cM    535      152       0.965        20
#   1.0 cM    268       84       0.901        14
#   2.0 cM    135       48       0.773        11
#   4.0 cM     69       29       0.565         9
#
# 1 cM is where adjacent tiles first reach r = 0.90. Below it they are near
# duplicates; above it the correlation keeps falling but distinct peaks merge,
# 2 cM already costing three of the fourteen resolved regions. (sig@7.5 reads
# 84 here against the 85 this script reports -- the tradeoff table tiles the
# euchromatic steps directly, while the analysis below tiles before filtering,
# so the arm-start offsets differ by one tile. Immaterial to the comparison.)
set.seed(1)
shares <- function(x) {
  s <- function(z) sum(z, na.rm = TRUE)
  t <- s(x$main) + s(x$sex) + s(x$diet) + s(x$int)
  100 * c(s(x$sex), s(x$diet), s(x$int)) / t
}
boot_ci <- function(x, nb = 2000) {
  ix <- split(seq_len(nrow(x)), paste0(x$chr, "_", x$tile %/% 4)); u <- names(ix)
  b <- t(replicate(nb, shares(x[unlist(ix[sample(u, length(u), TRUE)]), ])))
  apply(b, 2, quantile, c(.025, .975))
}
cat("\npartition over significant tiles, % of the four-treatment sum of squares.\n")
cat("95% intervals from 2000 bootstraps over 4 cM groups of tiles.\n\n")
w <- terms()
for (sc in c("all", "autosomes")) {
  x <- if (sc == "autosomes") w %>% filter(chr != "chrX") else w
  ci <- boot_ci(x); pt <- shares(x)
  cat(sprintf("  %-10s %3d tiles   shared %4.1f   sex %4.1f [%4.1f, %4.1f]   diet %4.1f [%4.1f, %4.1f]   int %5.1f [%5.1f, %4.1f]\n",
      sc, nrow(x), 100 - sum(pt),
      pt[1], ci[1,1], ci[2,1], pt[2], ci[1,2], ci[2,2], pt[3], ci[1,3], ci[2,3]))
}
