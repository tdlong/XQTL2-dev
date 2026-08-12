# varcomp_H2.R — partition Cutler H2 at each position.
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

TERMS     <- c("sex", "diet", "sex:diet")
TERM_COLS <- c(sex = "#1F78B4", diet = "#33A02C", `sex:diet` = "#E31A1C")
chr_order  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
chr_labels <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")

if (!file.exists(IN)) stop("input not found: ", IN)
long <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble()

w <- long %>%
  select(chr, pos, sex, sugar, half, Cutl_H2) %>%
  unite("cell", sugar, sex) %>%
  pivot_wider(names_from = c(cell, half), values_from = Cutl_H2, names_glue = "{cell}_{half}") %>%
  drop_na()

vc <- w %>% transmute(
  chr, pos,
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
         pct_main = 100 * main / tot, pct_sex = 100 * sex / tot,
         pct_diet = 100 * diet / tot, pct_int = 100 * `sex:diet` / tot,
         H2 = ybar)

write.table(vc %>% select(chr, pos, H2, ss_main, ss_sex, ss_diet, ss_int, ss_rep,
                          ms_rep, main, sex, diet, `sex:diet`,
                          pct_main, pct_sex, pct_diet, pct_int),
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
BIN <- 1e5
plot_df <- vc %>%
  mutate(chr = factor(chr, levels = chr_order), bin = (pos %/% BIN) * BIN) %>%
  filter(!is.na(chr)) %>%
  group_by(chr, bin) %>%
  summarise(sex = mean(pct_sex), diet = mean(pct_diet),
            `sex:diet` = mean(pct_int), H2 = mean(H2), .groups = "drop") %>%
  pivot_longer(all_of(TERMS), names_to = "term", values_to = "pct") %>%
  mutate(term = factor(term, levels = TERMS), pos_mb = bin / 1e6)

chr_lab_df <- plot_df %>% distinct(chr) %>% mutate(lab = chr_labels[as.character(chr)])

p <- ggplot(plot_df, aes(pos_mb, pct, colour = term)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.45) +
  facet_wrap(~ chr, ncol = 1, scales = "free_x") +
  geom_text(data = chr_lab_df, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.3, vjust = 1.4, size = 3, colour = "grey20",
            fontface = "bold", inherit.aes = FALSE) +
  scale_colour_manual(values = TERM_COLS, name = NULL) +
  scale_x_continuous(expand = expansion(0), minor_breaks = seq(0, 40, 1)) +
  labs(x = "Position (Mb)", y = "% of the H2 partition at that position",
       title = "How much of the heritability at each position is attributable to sex and to diet",
       subtitle = paste("Percent of n*ybar^2 + sex + diet + sex:diet, each corrected by the replicate error.",
                        "The remainder is the main term.\nThe main term is inflated -- the H2 floor sits in it",
                        "-- but these three are contrasts between treatments, so the floor cancels.")) +
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

cat("\nWhere sex accounts for the most:\n")
vc %>% slice_max(pct_sex, n = 5) %>%
  transmute(chr, pos_mb = round(pos/1e6, 2), H2 = round(H2, 2),
            `% sex` = round(pct_sex, 2), `% diet` = round(pct_diet, 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nWhere diet accounts for the most:\n")
vc %>% slice_max(pct_diet, n = 5) %>%
  transmute(chr, pos_mb = round(pos/1e6, 2), H2 = round(H2, 2),
            `% sex` = round(pct_sex, 2), `% diet` = round(pct_diet, 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)
