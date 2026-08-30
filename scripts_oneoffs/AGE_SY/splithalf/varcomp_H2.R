# varcomp_H2.R — partition Cutler h2 at each position into main, sex and diet.
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript scripts_oneoffs/AGE_SY/splithalf/varcomp_H2.R [long.txt.gz]
#
# Input is the long file from gather_splithalf_H2.R: 8 measures per window,
# 4 treatments x 2 halves (odd and even replicates).
#
# UNITS. h2 is a PERCENTAGE -- Heritability() multiplies by 200 = 2 x 100 -- so
# a reported 0.7 is 0.7% of the variance at that window. Every term below is in
# percentage points.
#
# MODEL. A balanced 2x2 with 2 replicates per cell, decomposed uncorrected:
#
#   sum(y^2) = n*ybar^2 + SS_sex + SS_diet + SS_sex:diet + SS_rep
#
# n*ybar^2 = 8*ybar^2 and the treatment terms carry 1 df each; SS_rep, from the
# within-cell odd-vs-even differences, carries 4 and is pure error. Because the
# df differ, the error is subtracted on MEAN squares: each term is MS - SS_rep/4.
# Checked against aov() at chr3L 9,400,000; the parts sum to sum(y^2) to 9e-14.
#
# Terms are also given in h2 units as sign * sqrt(|MS - MS_rep| / 8), since every
# term satisfies SS = 8 * deviation^2. They are NOT floored at zero: a term below
# the replicate error goes negative, and that is the signal that nothing is there.
#
# THE h2 FLOOR. h2 squares an estimated quantity, so E[h2] = true + b, and b does
# not shrink with replication -- it is a within-replicate bias identical in every
# replicate. XQTL2 #34 reports b per window. The reported b is well calibrated
# where it is small but badly overstated where it is large, which drives corrected
# h2 to -6 at a few Mb. It is therefore recalibrated against windows the Wald test
# calls null: where the frequencies did not move, true h2 ~ 0, so the observed h2
# IS the floor. An isotonic fit of h2 on b over those windows is applied
# everywhere. The fitted curve is nearly flat, so in practice this subtracts a
# near-constant ~0.55 and the spatial structure in the reported b is not real.

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
IN   <- if (length(args) >= 1) args[1] else
  "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
# Second argument sets the output, so a variant input (e.g. the no-8-9 halves)
# does not overwrite the 12-replicate partition.
OUT  <- if (length(args) >= 2) args[2] else
        "process/AGE_SY_splithalf/H2_varcomp_by_window.txt.gz"

WALD_NULL <- 2      # Wald below this: frequencies did not move, so true h2 ~ 0
POOL_BP   <- 0      # >0 pools MS_rep over this span before subtracting

TERMS <- c("main", "sex", "diet", "sex:diet")

if (!file.exists(IN)) stop("input not found: ", IN)
long <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble()

# ── floor: calibrate the reported bias against Wald-null windows ─────────────
# No floor calibration. H2 from hap_scan (XQTL2 #40) averages replicates before
# squaring and takes its correction from them, so it carries no squaring bias to
# calibrate away -- the isotonic fit of h2 on b that used to live here was
# repairing an estimator that no longer exists. The split-half error term below
# is a different thing and stays: it removes the VARIANCE of the h2 estimate,
# which the within-window correction says nothing about.
# H2 from hap_scan already carries its own bias correction, measured from the
# replicates, so nothing is subtracted here. The value entering the decomposition
# is the heritability as reported.
long <- long %>% mutate(H2v = H2)

# ── one row per window, four cell means and the within-cell error ────────────
wald_max <- long %>%
  group_by(chr, pos, sex, sugar) %>% summarise(w = mean(Wald_log10p), .groups = "drop") %>%
  group_by(chr, pos) %>% summarise(wald_max = max(w), .groups = "drop")

vc <- long %>%
  select(chr, pos, sex, sugar, half, H2v) %>%
  unite("cell", sugar, sex) %>%
  pivot_wider(names_from = c(cell, half), values_from = H2v, names_glue = "{cell}_{half}") %>%
  drop_na() %>%
  left_join(wald_max, by = c("chr", "pos")) %>%
  transmute(
    chr, pos, wald_max,
    a  = (SY10_F_odd + SY10_F_even) / 2,     # SY10 female
    b  = (SY20_F_odd + SY20_F_even) / 2,     # SY20 female
    cc = (SY10_M_odd + SY10_M_even) / 2,     # SY10 male
    d  = (SY20_M_odd + SY20_M_even) / 2,     # SY20 male
    ss_rep = ((SY10_F_odd - SY10_F_even)^2 + (SY20_F_odd - SY20_F_even)^2 +
              (SY10_M_odd - SY10_M_even)^2 + (SY20_M_odd - SY20_M_even)^2) / 2) %>%
  mutate(
    ybar = (a + b + cc + d) / 4,
    F_ = (a + b) / 2, M_ = (cc + d) / 2, S10 = (a + cc) / 2, S20 = (b + d) / 2,
    ss_main = 8 * ybar^2,
    ss_sex  = 4 * ((F_  - ybar)^2 + (M_  - ybar)^2),
    ss_diet = 4 * ((S10 - ybar)^2 + (S20 - ybar)^2),
    ss_int  = 2 * ((a  - F_ - S10 + ybar)^2 + (b - F_ - S20 + ybar)^2 +
                   (cc - M_ - S10 + ybar)^2 + (d - M_ - S20 + ybar)^2),
    ms_rep  = ss_rep / 4)

if (POOL_BP > 0)
  vc <- vc %>% group_by(region = paste0(chr, "_", pos %/% POOL_BP)) %>%
    mutate(ms_rep = mean(ms_rep)) %>% ungroup() %>% select(-region)

# ── subtract the pure error; convert to h2 units ─────────────────────────────
vc <- vc %>%
  mutate(main = ss_main - ms_rep, sex = ss_sex - ms_rep,
         diet = ss_diet - ms_rep, `sex:diet` = ss_int - ms_rep,
         # SS = 8 * deviation^2, so sqrt(SS/8) returns h2 units. main is signed
         # by ybar, not by the term: 8*ybar^2 is a square, so an over-corrected
         # window with strongly negative h2 would otherwise read as a large
         # positive shared effect.
         mainH2 = sign(ybar)       * sqrt(abs(main) / 8),
         sexH2  = sign(sex)        * sqrt(abs(sex) / 8),
         dietH2 = sign(diet)       * sqrt(abs(diet) / 8),
         intH2  = sign(`sex:diet`) * sqrt(abs(`sex:diet`) / 8),
         dir_sex  = sign(M_ - F_),      # +1 male higher
         dir_diet = sign(S20 - S10),    # +1 SY20 higher
         H2 = ybar)

write.table(vc %>% select(chr, pos, wald_max, H2, ms_rep,
                          all_of(TERMS), mainH2, sexH2, dietH2, intH2,
                          dir_sex, dir_diet),
            gzfile(OUT), row.names = FALSE, quote = FALSE, sep = "\t")
cat("Wrote:", OUT, "\n")

# ── report ───────────────────────────────────────────────────────────────────
partition <- function(d, label) {
  s <- d %>% summarise(across(all_of(TERMS), sum))
  t <- sum(unlist(s))
  cat(sprintf("%-26s main %5.1f%%  sex %5.1f%%  diet %5.1f%%  sex:diet %5.1f%%  (n=%d)\n",
              label, 100*s$main/t, 100*s$sex/t, 100*s$diet/t, 100*s$`sex:diet`/t, nrow(d)))
}
cat("\nVariance partition of Cutler h2\n")
partition(vc %>% filter(chr != "chrX"), "autosomes (headline)")
partition(vc,                            "whole genome incl. X")
partition(vc %>% filter(chr == "chrX"),  "X only")
partition(vc %>% filter(chr == "chr3L", pos > 9.3e6, pos < 9.6e6), "chr3L QTL 9.3-9.6 Mb")

cat("\nBelow the replicate error at:\n")
for (tm in TERMS) cat(sprintf("  %-9s %5.1f%% of windows\n", tm, 100*mean(vc[[tm]] < 0)))

cat("\nLargest terms (h2 percentage points):\n")
for (v in c("sexH2", "dietH2")) {
  cat(" ", v, "\n")
  vc %>% slice_max(.data[[v]], n = 3) %>%
    transmute(chr, Mb = round(pos/1e6, 2), H2 = round(H2, 2), wald = round(wald_max, 1),
              mainH2 = round(mainH2, 2), sexH2 = round(sexH2, 3), dietH2 = round(dietH2, 3)) %>%
    as.data.frame() %>% print(row.names = FALSE)
}
