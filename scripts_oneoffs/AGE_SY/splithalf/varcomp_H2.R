# varcomp_H2.R — partition Cutler H2 at each position.
#
# UNITS: H2 is a PERCENTAGE. Heritability() multiplies by 200 (= 2 x 100), so a
# reported 0.7 means 0.7% of the variance at that window, and the reported bias
# carries the same factor. For a polygenic trait that is the expected scale per
# 75 kb window; the chr3L peak at 4-5.5% is the outlier. All the terms below are
# therefore in percentage points, not proportions.
#
# H2 is the response. Eight measures per position: 4 treatments x 2 halves, the
# halves being pure replicates. The uncorrected sum of squares decomposes as
#
#   sum(y^2) = n*ybar^2 + SS_sex + SS_diet + SS_sex:diet + SS_rep
#
# with, for this balanced design (a = b = 2 levels, n = 2 replicates per cell):
#
#   n*ybar^2      = 8 * ybar^2                                    1 df
#   SS_sex        = 4 * sum_s (ybar_s - ybar)^2                    1 df
#   SS_diet       = 4 * sum_d (ybar_d - ybar)^2                    1 df
#   SS_sex:diet   = 2 * sum_sd (ybar_sd - ybar_s - ybar_d + ybar)^2  1 df
#   SS_rep        = sum_sdr (y_sdr - ybar_sd)^2                    4 df
#
# The replicate term is the pure error. Since it carries 4 df and the others 1,
# the subtraction is done on MEAN squares: MS_rep = SS_rep/4, and each term
# becomes MS - MS_rep. Percentages are of the sum of the four corrected terms.
#
# WHAT THE CORRECTION DOES AND DOES NOT DO
#   Writing each replicate as  H2_r = true + b + e_r,  MS_rep estimates Var(e_r)
#   and subtracting it removes exactly that. It does NOT remove b, the reason H2
#   is positive everywhere: b comes from Affect^2 squaring a noisy quantity
#   WITHIN each replicate, so every replicate carries the same b, the replicates
#   agree about it, and the error term is blind to it. b sits in the main term,
#   which is therefore inflated. Measuring b needs a null with no true
#   heritability -- controls split into two pseudo-groups.
#
#   sex, diet and sex:diet are contrasts BETWEEN cells, so b cancels exactly.
#   Those three are clean.
#
# POOL_BP > 0 pools MS_rep over a local run of windows before subtracting; 4 df
# per window is noisy, and windows step 5 kb under a 75 kb haplotype window so
# neighbours are near-duplicates. 0 (the default) uses each window's own MS_rep.
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript scripts_oneoffs/AGE_SY/splithalf/varcomp_H2.R [long.txt.gz]

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
IN   <- if (length(args) >= 1) args[1] else
  "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
OUT  <- "process/AGE_SY_splithalf/H2_varcomp_by_window.txt.gz"
FIG  <- "figures/AGE_SY_H2_varcomp.png"
POOL_BP <- 0

# Wald tests whether the haplotype frequencies changed at all. Where it is not
# significant, H2 is computed from differences we have already called noise --
# and because H2 squares those differences, noise still yields a positive H2.
# So the partition is only reported where at least one treatment shows a real
# frequency change. 6 is the threshold used on the 4-scan Wald figure.
WALD_MIN <- 6
#
# How precisely is H2 measured, as a function of Wald? Binning on ONE half's
# Wald and reading the other half's H2 (binning on the mean of both conditions
# on a collider and manufactures negative correlations), and counting
# INDEPENDENT REGIONS rather than windows -- windows step 5 kb under a 75 kb
# haplotype window, so window counts are inflated ~50x:
#
#   Wald    regions  distinct Mb   |odd-even|   rel err
#   <2        414       110          0.169       31%
#   4-6       124        80          0.355       42%
#   10-15      60        39          0.587       39%
#   15-25      23        23          0.864       55%
#   25-40       8         3          0.551       22%
#   40-60       6         2          0.246        6%
#   >60         3         2          0.057        1%
#
# Over the range with enough independent regions to support a statement
# (Wald < 25), relative error is 30-55% and does NOT improve with Wald;
# absolute error tracks H2, which is the multiplicative-error pattern. The
# apparent collapse above 25 is not a trend: every window there is chr3L 8-9 Mb,
# almost all of it male. So Wald > 6 marks where H2 MEANS something, not where
# it is measured well. In this experiment H2 is imprecise nearly everywhere
# except the chr3L peak.

TERMS     <- c("sex", "diet", "sex:diet")
TERM_COLS <- c(sex = "#1F78B4", diet = "#33A02C", `sex:diet` = "#E31A1C")
chr_order  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
chr_labels <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")

if (!file.exists(IN)) stop("input not found: ", IN)
long <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble()

# XQTL2 #34: H2 has an upward floor from squaring a noisy Affect, and the
# pipeline now reports it. Subtract it from EACH of the 8 measures before
# decomposing -- the floor lives in the main term, which is the denominator of
# every percentage below, so leaving it in makes sex and diet read low. It also
# differs between cells (it tracks depth and fly counts), so a per-cell
# subtraction cleans the treatment terms slightly too.
HAVE_BIAS <- "Cutl_H2_bias" %in% names(long) && !all(is.na(long$Cutl_H2_bias))
if (HAVE_BIAS) {
  cat(sprintf("Subtracting the reported H2 floor: median bias %.3f against median H2 %.3f\n",
              median(long$Cutl_H2_bias, na.rm = TRUE), median(long$Cutl_H2, na.rm = TRUE)))
  if ("Cutl_clamp_frac" %in% names(long))
    cat(sprintf("  median penetrance-clamp fraction %.3f (high values = bias less trustworthy)\n",
                median(long$Cutl_clamp_frac, na.rm = TRUE)))
  long <- long %>% mutate(H2v = Cutl_H2 - Cutl_H2_bias)
} else {
  cat("NOTE: no Cutl_H2_bias column -- these scans predate XQTL2 #34.\n")
  cat("      The main term keeps its floor, so the sex and diet percentages\n")
  cat("      below are LOWER BOUNDS. Re-run the scans to fix this.\n")
  long <- long %>% mutate(H2v = Cutl_H2)
}

# Strongest frequency change at each window, across the four treatments.
wald_max <- long %>%
  group_by(chr, pos, sex, sugar) %>%
  summarise(w = mean(Wald_log10p), .groups = "drop") %>%
  group_by(chr, pos) %>% summarise(wald_max = max(w), .groups = "drop")

w <- long %>%
  select(chr, pos, sex, sugar, half, H2v) %>%
  unite("cell", sugar, sex) %>%
  pivot_wider(names_from = c(cell, half), values_from = H2v, names_glue = "{cell}_{half}") %>%
  drop_na() %>%
  left_join(wald_max, by = c("chr", "pos"))

cat(sprintf("Windows with Wald > %g in at least one treatment: %d of %d (%.1f%%)\n",
            WALD_MIN, sum(w$wald_max > WALD_MIN), nrow(w),
            100 * mean(w$wald_max > WALD_MIN)))

vc <- w %>% transmute(
  chr, pos, wald_max,
  a  = (SY10_F_odd + SY10_F_even) / 2,      # cell means
  b  = (SY20_F_odd + SY20_F_even) / 2,
  cc = (SY10_M_odd + SY10_M_even) / 2,
  d  = (SY20_M_odd + SY20_M_even) / 2,
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
    ms_rep_w = ss_rep / 4)

vc <- if (POOL_BP > 0) {
  vc %>% mutate(region = paste0(chr, "_", pos %/% POOL_BP)) %>%
    group_by(region) %>% mutate(ms_rep = mean(ms_rep_w)) %>% ungroup()
} else vc %>% mutate(ms_rep = ms_rep_w)

# MS = SS for the 1 df terms; subtract the pure error from each
vc <- vc %>%
  mutate(main = ss_main - ms_rep, sex = ss_sex - ms_rep,
         diet = ss_diet - ms_rep, `sex:diet` = ss_int - ms_rep,
         tot  = main + sex + diet + `sex:diet`,
         # Back to H2 UNITS (percentage points). Every term satisfies
         # SS = 8*deviation^2, so the deviation is sqrt(SS/8) -- comparable to the
         # four treatment curves, which sums of squares are not.
         #
         # NOT floored at zero. A term goes negative when it sits below the
         # replicate error, and that is the diagnostic worth seeing: it says the
         # correction is larger than the signal, i.e. there is nothing here. So
         # the signed root keeps the sign of the COMPONENT. Direction of the
         # contrast (which sex, which diet is higher) is a separate column.
         sexH2  = sign(sex)        * sqrt(abs(sex) / 8),
         dietH2 = sign(diet)       * sqrt(abs(diet) / 8),
         intH2  = sign(`sex:diet`) * sqrt(abs(`sex:diet`) / 8),
         mainH2 = sign(main)       * sqrt(abs(main) / 8),
         dir_sex  = sign(M_ - F_),      # +1 male higher
         dir_diet = sign(S20 - S10),    # +1 SY20 higher
         pct_main = 100 * main / tot, pct_sex = 100 * sex / tot,
         pct_diet = 100 * diet / tot, pct_int = 100 * `sex:diet` / tot,
         H2 = ybar)

write.table(vc %>% select(chr, pos, wald_max, H2, ss_main, ss_sex, ss_diet, ss_int,
                          ss_rep, ms_rep, main, sex, diet, `sex:diet`,
                          pct_main, pct_sex, pct_diet, pct_int,
                          mainH2, sexH2, dietH2, intH2, dir_sex, dir_diet),
            gzfile(OUT), row.names = FALSE, quote = FALSE, sep = "\t")
cat("Wrote:", OUT, "\n")
cat(if (POOL_BP > 0) sprintf("MS_rep pooled over %.0f kb\n", POOL_BP/1e3)
    else "MS_rep taken per window (4 df)\n")

report <- function(df, label) {
  s <- df %>% summarise(across(c(H2, ss_main, ss_sex, ss_diet, ss_int, ss_rep,
                                 ms_rep, main, sex, diet, `sex:diet`), mean))
  cat("\n", label, "   H2 = ", round(s$H2, 3), "   MS_rep = ", round(s$ms_rep, 4), "\n", sep = "")
  tibble(term = c("main", "sex", "diet", "sex:diet", "rep"),
         df   = c(1, 1, 1, 1, 4),
         SS   = c(s$ss_main, s$ss_sex, s$ss_diet, s$ss_int, s$ss_rep)) %>%
    mutate(MS = SS / df,
           `MS-MSrep` = ifelse(term == "rep", NA, MS - s$ms_rep),
           `%` = ifelse(term == "rep", NA,
                        100 * `MS-MSrep` / (s$main + s$sex + s$diet + s$`sex:diet`)),
           across(where(is.numeric), ~ round(.x, 4))) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

report(vc %>% filter(chr == "chr3L", pos == 9355000), "chr3L 9,355,000 (peak)")

# ── the picture: how big are sex and diet relative to the whole, by position ──
# Plot the ABSOLUTE terms, whole genome, no mask. Each is MS - MS_rep, so where
# nothing is happening it sits at zero by construction -- the subtraction does
# the masking. Percentages cannot be plotted this way: a percentage of a tiny
# total is still a percentage, which is what forced a Wald cut before.
BIN <- 1e5
plot_df <- vc %>%
  mutate(chr = factor(chr, levels = chr_order), bin = (pos %/% BIN) * BIN) %>%
  filter(!is.na(chr)) %>%
  group_by(chr, bin) %>%
  summarise(sex = mean(sexH2), diet = mean(dietH2), `sex:diet` = mean(intH2),
            H2 = mean(H2), .groups = "drop") %>%
  pivot_longer(all_of(TERMS), names_to = "term", values_to = "pct") %>%
  mutate(term = factor(term, levels = TERMS), pos_mb = bin / 1e6)


chr_lab_df <- plot_df %>% distinct(chr) %>% mutate(lab = chr_labels[as.character(chr)])

p <- ggplot(plot_df, aes(pos_mb, pct, colour = term)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ chr, ncol = 1, scales = "free_x") +
  geom_text(data = chr_lab_df, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.3, vjust = 1.4, size = 3, colour = "grey20",
            fontface = "bold", inherit.aes = FALSE) +
  scale_colour_manual(values = TERM_COLS, name = NULL) +
  scale_x_continuous(expand = expansion(0), minor_breaks = seq(0, 40, 1)) +
  labs(x = "Position (Mb)",
       y = expression("Cutler" ~ H^2 ~ "(percentage points): deviation, replicate error subtracted"),
       title = if (HAVE_BIAS)
           "Variance in H2 attributable to sex and to diet, whole genome"
         else
           "PRELIMINARY - H2 floor NOT subtracted, do not interpret (see XQTL2 #34)",
       subtitle = paste0("Signed sqrt((MS - MS_rep)/8): back in H2 percentage points, comparable to the four treatment curves.\n",
                         "NEGATIVE means the term is below the replicate error -- nothing there. Not floored at zero.")) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top",
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        strip.background = element_blank(), strip.text = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing = unit(0.4, "lines"))

dir.create("figures", showWarnings = FALSE)
png(FIG, width = 8.5, height = 6.5, units = "in", res = 150)
print(p)
invisible(dev.off())
cat("\nWrote:", FIG, "\n")

# Rankings in H2 UNITS, not percentages: once the floor is removed the total is
# near zero over most of the genome, so a percentage of it is meaningless
# (23,000% and worse). Percentages are only sane where H2 is well above zero.
cat("\nLargest sex terms (H2 percentage points):\n")
vc %>% slice_max(sexH2, n = 5) %>%
  transmute(chr, Mb = round(pos/1e6, 2), H2 = round(H2, 2), wald = round(wald_max, 1),
            mainH2 = round(mainH2, 2), sexH2 = round(sexH2, 3), dietH2 = round(dietH2, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nLargest diet terms (H2 percentage points):\n")
vc %>% slice_max(dietH2, n = 5) %>%
  transmute(chr, Mb = round(pos/1e6, 2), H2 = round(H2, 2), wald = round(wald_max, 1),
            mainH2 = round(mainH2, 2), sexH2 = round(sexH2, 3), dietH2 = round(dietH2, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
