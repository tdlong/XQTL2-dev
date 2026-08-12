# varcomp_H2.R — H2 is the response. Partition its variance at each position.
#
# Eight measures per position: 4 treatments x 2 halves. The halves are pure
# replicates, so at each window
#
#   Source      df   E[MS]
#   sex          1   s2_e + 2*s2_int + 4*s2_sex
#   diet         1   s2_e + 2*s2_int + 4*s2_diet
#   sex:diet     1   s2_e + 2*s2_int
#   rep(cell)    4   s2_e
#
# solved for the components:
#
#   s2_int  = (MS_int  - MS_e)   / 2
#   s2_sex  = (MS_sex  - MS_int) / 4
#   s2_diet = (MS_diet - MS_int) / 4
#   s2_e    =  MS_e
#
# Two views are reported. The sums of squares decompose the total exactly with
# no assumptions (SS_total = SS_sex + SS_diet + SS_int + SS_error), and the
# variance components above additionally strip the replicate noise out of the
# treatment terms.
#
# MS_e on 4 df is noisy per window, and windows step 5 kb under a 75 kb
# haplotype window so neighbours are near-duplicates; it is pooled over POOL_BP.
# That pools the ERROR only.
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript scripts_oneoffs/AGE_SY/splithalf/varcomp_H2.R [long.txt.gz]

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
IN   <- if (length(args) >= 1) args[1] else
  "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
OUT  <- "process/AGE_SY_splithalf/H2_varcomp_by_window.txt.gz"
FIG  <- "figures/AGE_SY_H2_varcomp.png"
POOL_BP <- 5e5

SRC      <- c("sex", "diet", "sex:diet", "error")
SRC_COLS <- c(sex = "#1F78B4", diet = "#33A02C", `sex:diet` = "#E31A1C", error = "grey70")
chr_order  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
chr_labels <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")

if (!file.exists(IN)) stop("input not found: ", IN)
long <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble()

w <- long %>%
  select(chr, pos, sex, sugar, half, Cutl_H2) %>%
  unite("cell", sugar, sex) %>%
  pivot_wider(names_from = c(cell, half), values_from = Cutl_H2, names_glue = "{cell}_{half}") %>%
  drop_na()

av <- w %>% transmute(
  chr, pos,
  a  = (SY10_F_odd + SY10_F_even) / 2,
  b  = (SY20_F_odd + SY20_F_even) / 2,
  cc = (SY10_M_odd + SY10_M_even) / 2,
  d  = (SY20_M_odd + SY20_M_even) / 2,
  ss_error = ((SY10_F_odd - SY10_F_even)^2 + (SY20_F_odd - SY20_F_even)^2 +
              (SY10_M_odd - SY10_M_even)^2 + (SY20_M_odd - SY20_M_even)^2) / 2) %>%
  mutate(region = paste0(chr, "_", pos %/% POOL_BP)) %>%
  group_by(region) %>% mutate(ss_error = mean(ss_error), n_pool = n()) %>% ungroup()

vc <- av %>%
  mutate(
    grand = (a + b + cc + d) / 4,
    F_ = (a + b) / 2, M_ = (cc + d) / 2,
    S10 = (a + cc) / 2, S20 = (b + d) / 2,
    ss_sex  = 2 * 2 * ((F_ - grand)^2 + (M_ - grand)^2),
    ss_diet = 2 * 2 * ((S10 - grand)^2 + (S20 - grand)^2),
    ss_int  = 2 * ((a - F_ - S10 + grand)^2 + (b - F_ - S20 + grand)^2 +
                   (cc - M_ - S10 + grand)^2 + (d - M_ - S20 + grand)^2),
    ss_total = ss_sex + ss_diet + ss_int + ss_error,
    ms_sex = ss_sex, ms_diet = ss_diet, ms_int = ss_int,   # 1 df each
    ms_e   = ss_error / 4,
    vc_int  = pmax((ms_int  - ms_e)   / 2, 0),
    vc_sex  = pmax((ms_sex  - ms_int) / 4, 0),
    vc_diet = pmax((ms_diet - ms_int) / 4, 0),
    vc_err  = ms_e,
    vc_tot  = vc_sex + vc_diet + vc_int + vc_err,
    H2 = grand)

write.table(vc %>% select(chr, pos, H2, ss_sex, ss_diet, ss_int, ss_error,
                          vc_sex, vc_diet, vc_int, vc_err),
            gzfile(OUT), row.names = FALSE, quote = FALSE, sep = "\t")
cat("Wrote:", OUT, "\n")
cat(sprintf("MS_error pooled over %.0f kb (%.0f windows, 4 df each)\n\n",
            POOL_BP / 1e3, median(vc$n_pool)))

report <- function(df, label) {
  s <- df %>% summarise(across(c(H2, starts_with("ss_"), starts_with("vc_")), mean))
  cat(label, "\n", sep = "")
  cat(sprintf("  H2 here (mean of the 8) = %.3f\n", s$H2))
  cat("             sum of squares        variance component\n")
  ssv <- c(s$ss_sex, s$ss_diet, s$ss_int, s$ss_error)
  vcv <- c(s$vc_sex, s$vc_diet, s$vc_int, s$vc_err)
  for (i in seq_along(SRC))
    cat(sprintf("  %-9s %8.3f  %5.1f%%      %8.4f  %5.1f%%\n",
                SRC[i], ssv[i], 100*ssv[i]/sum(ssv), vcv[i], 100*vcv[i]/sum(vcv)))
  cat(sprintf("  %-9s %8.3f            %8.4f\n\n", "total", sum(ssv), sum(vcv)))
}

report(vc %>% filter(chr == "chr3L", pos > 9.3e6, pos < 9.5e6), "chr3L 9.3-9.5 Mb")
report(vc, "genome-wide average")

# ── F tests, so "sex does not matter here" can be said rather than eyeballed ──
vc <- vc %>%
  mutate(F_sex  = ss_sex  / ms_e, F_diet = ss_diet / ms_e, F_int = ss_int / ms_e,
         p_sex  = pf(F_sex,  1, 4, lower.tail = FALSE),
         p_diet = pf(F_diet, 1, 4, lower.tail = FALSE),
         p_int  = pf(F_int,  1, 4, lower.tail = FALSE))

BIN <- 2e5
stacked <- vc %>%
  mutate(chr = factor(chr, levels = chr_order), bin = (pos %/% BIN) * BIN) %>%
  filter(!is.na(chr)) %>%
  group_by(chr, bin) %>%
  summarise(sex = mean(vc_sex), diet = mean(vc_diet),
            `sex:diet` = mean(vc_int), error = mean(vc_err), .groups = "drop") %>%
  pivot_longer(all_of(SRC), names_to = "source", values_to = "var") %>%
  mutate(source = factor(source, levels = SRC), pos_mb = bin / 1e6)

chr_lab_df <- stacked %>% distinct(chr) %>% mutate(lab = chr_labels[as.character(chr)])

# The variance stack alone cannot separate "strong QTL, same in all four
# treatments" from "nothing here" -- both have near-zero variance among the 8
# measures. So the H2 LEVEL is the line, and the rug underneath marks where sex
# or diet is significant. High line + no rug = heritability that sex and diet
# do not touch.
HET <- tribble(~chr, ~eu_start, ~eu_end,
  "chrX", 2.5, 21.2, "chr2L", 0.5, 22.9, "chr2R", 1.3, 25.1,
  "chr3L", 0.7, 24.0, "chr3R", 4.5, 32.0)

lvl <- vc %>%
  mutate(chr = factor(chr, levels = chr_order), pos_mb = pos / 1e6,
         bin = (pos %/% BIN) * BIN) %>%
  filter(!is.na(chr)) %>%
  group_by(chr, bin) %>%
  summarise(H2 = mean(H2), pos_mb = mean(pos_mb),
            sex_sig = mean(p_sex < 0.05), diet_sig = mean(p_diet < 0.05),
            .groups = "drop")

bg <- median(vc$H2)
ymax <- max(lvl$H2) * 1.08
rug <- bind_rows(
  lvl %>% filter(sex_sig  > 0.5) %>% transmute(chr, pos_mb, term = "sex matters here",
                                               ymin = -0.28*bg, ymax = -0.06*bg),
  lvl %>% filter(diet_sig > 0.5) %>% transmute(chr, pos_mb, term = "diet matters here",
                                               ymin = -0.52*bg, ymax = -0.30*bg))

het_rects <- HET %>% mutate(chr = factor(chr, levels = chr_order)) %>%
  left_join(lvl %>% group_by(chr) %>% summarise(xmax = max(pos_mb), .groups="drop"), by="chr") %>%
  { bind_rows(transmute(., chr, xmin = 0, xmax2 = eu_start),
              transmute(., chr, xmin = eu_end, xmax2 = xmax)) } %>% rename(xmax = xmax2)

p <- ggplot(lvl, aes(pos_mb, H2)) +
  geom_rect(data = het_rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey88", alpha = 0.6, colour = NA, inherit.aes = FALSE) +
  geom_hline(yintercept = bg, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
  geom_area(fill = "grey75", colour = NA) +
  geom_line(linewidth = 0.4) +
  geom_rect(data = rug, aes(xmin = pos_mb - BIN/2e6, xmax = pos_mb + BIN/2e6,
                            ymin = ymin, ymax = ymax, fill = term),
            inherit.aes = FALSE) +
  facet_wrap(~ chr, ncol = 1, scales = "free_x") +
  geom_text(data = chr_lab_df, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.3, vjust = 1.4, size = 3, colour = "grey20",
            fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = c("sex matters here" = "#1F78B4",
                               "diet matters here" = "#33A02C"), name = NULL) +
  scale_x_continuous(expand = expansion(0), minor_breaks = seq(0, 40, 1)) +
  scale_y_continuous(limits = c(-0.55*bg, ymax), expand = expansion(0)) +
  labs(x = "Position (Mb)", y = expression("Cutler" ~ H^2),
       title = "Where is the heritability, and where do sex or diet actually account for it",
       subtitle = paste0("Dashed line is the genome-wide background (", round(bg, 2),
                         "). Grey shading is heterochromatin.\n",
                         "Bars below the axis mark 1 Mb bins where sex (blue) or diet (green) is ",
                         "significant, F on 1 and 4 df, p<0.05 in most windows.\n",
                         "A tall curve with no bar under it is heritability that sex and diet do not account for.")) +
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
cat("Wrote:", FIG, "\n")

cat("\n1 Mb bins, H2 level vs whether sex/diet are significant:\n")
lvl %>%
  mutate(level = ifelse(H2 > quantile(lvl$H2, 0.75), "H2 elevated", "H2 at background"),
         terms = case_when(sex_sig > 0.5 & diet_sig > 0.5 ~ "sex and diet",
                           sex_sig > 0.5 ~ "sex only", diet_sig > 0.5 ~ "diet only",
                           TRUE ~ "neither")) %>%
  count(level, terms) %>% pivot_wider(names_from = terms, values_from = n, values_fill = 0) %>%
  as.data.frame() %>% print(row.names = FALSE)
