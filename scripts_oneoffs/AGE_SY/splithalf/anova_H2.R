# anova_H2.R — the nested partition at each position.
#
#   sex + diet + sex:diet + rep %in% (sex:diet)
#
# The odd and even halves ARE the replicates: two independent estimates of H2
# per cell per window, so at every position there are 8 numbers, and the four
# within-cell odd-vs-even differences give a 4 df error term. No smoothing, no
# genome-wide covariance -- the noise is measured at the position itself.
#
# Effects are in DEVIATION form, so they add back as
#
#   H2(sugar, sex) = mean  +/- sex  +/- diet  +/- interaction
#
# i.e. where the four treatments read 2.8, 2.8, 4.9, 5.5 the effects read
# mean 4.0, sex 1.2, diet 0.15, interaction 0.15.
#
# Each effect is a contrast in the four cell means, each mean over n=2 halves,
# with coefficients +-1/4, so
#
#   var(effect) = MS_error * sum(coef^2) / n = MS_error / 8
#
# and the noise-corrected squared effect is  effect^2 - MS_error/8, which is
# what "how much of the variation here is due to sex" actually means. It can go
# negative, which is the honest reading of "nothing beyond noise".
#
# MS_error on 4 df is very noisy per window. Because windows step 5 kb under a
# 75 kb haplotype window, neighbours are near-duplicates, so MS_error is pooled
# over a local run of windows (POOL_BP) -- that is averaging the error term over
# a region, not smoothing the effects.
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript scripts_oneoffs/AGE_SY/splithalf/anova_H2.R [long_file.txt.gz]

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
IN   <- if (length(args) >= 1) args[1] else
  "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
OUT  <- "process/AGE_SY_splithalf/H2_anova_by_window.txt.gz"
FIG  <- "figures/AGE_SY_H2_anova.png"
POOL_BP <- 5e5     # region over which MS_error is pooled

TERMS     <- c("mean", "sex", "diet", "sex:diet")
TERM_COLS <- c(mean = "black", sex = "#1F78B4", diet = "#33A02C", `sex:diet` = "#E31A1C")
chr_order  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
chr_labels <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")

if (!file.exists(IN)) stop("input not found: ", IN)
long <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble()

w <- long %>%
  select(chr, pos, sex, sugar, half, Cutl_H2) %>%
  unite("cell", sugar, sex) %>%
  pivot_wider(names_from = c(cell, half), values_from = Cutl_H2,
              names_glue = "{cell}_{half}") %>%
  drop_na()

# ── cell means (over the two halves) and the within-cell error ───────────────
cm <- w %>% transmute(
  chr, pos,
  SY10_F = (SY10_F_odd + SY10_F_even) / 2,
  SY20_F = (SY20_F_odd + SY20_F_even) / 2,
  SY10_M = (SY10_M_odd + SY10_M_even) / 2,
  SY20_M = (SY20_M_odd + SY20_M_even) / 2,
  # SS_error = sum over the 4 cells of (odd-even)^2 / 2, on 4 df
  ms_error = ((SY10_F_odd - SY10_F_even)^2 + (SY20_F_odd - SY20_F_even)^2 +
              (SY10_M_odd - SY10_M_even)^2 + (SY20_M_odd - SY20_M_even)^2) / 2 / 4)

# Pool MS_error locally: 4 df per window is far too noisy, and neighbouring
# windows are near-duplicates anyway (5 kb steps under a 75 kb window).
cm <- cm %>%
  mutate(region = paste0(chr, "_", pos %/% POOL_BP)) %>%
  group_by(region) %>%
  mutate(ms_error_pooled = mean(ms_error), n_pool = n()) %>%
  ungroup()

eff <- cm %>%
  transmute(
    chr, pos, ms_error = ms_error_pooled,
    mean       = (SY10_F + SY20_F + SY10_M + SY20_M) / 4,
    sex        = (SY10_F + SY20_F - SY10_M - SY20_M) / 4,
    diet       = (SY20_F + SY20_M - SY10_F - SY10_M) / 4,
    `sex:diet` = (SY10_F - SY20_F - SY10_M + SY20_M) / 4)

# var(effect) = MS_error * sum(coef^2) / n  with coef = +-1/4, n = 2  ->  /8
eff <- eff %>%
  mutate(across(all_of(TERMS), ~ .x^2 - ms_error / 8, .names = "var_{.col}"))

write.table(eff, gzfile(OUT), row.names = FALSE, quote = FALSE, sep = "\t")
cat("Wrote:", OUT, "\n")

cat(sprintf("\nMS_error pooled over %.0f kb (%.0f windows per region, 4 df each)\n",
            POOL_BP / 1e3, median(cm$n_pool)))

cat("\nGenome-wide medians (H2 units):\n")
eff %>% summarise(across(all_of(TERMS), median)) %>%
  pivot_longer(everything(), names_to = "term", values_to = "effect") %>%
  left_join(eff %>% summarise(across(starts_with("var_"), median)) %>%
              pivot_longer(everything(), names_to = "term", values_to = "var_corrected") %>%
              mutate(term = sub("^var_", "", term)), by = "term") %>%
  mutate(term = factor(term, levels = TERMS), across(where(is.numeric), ~ round(.x, 4))) %>%
  arrange(term) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\nAt the chr3L peak (9.3-9.5 Mb):\n")
eff %>% filter(chr == "chr3L", pos > 9.3e6, pos < 9.5e6) %>%
  summarise(across(c(all_of(TERMS), ms_error), mean)) %>%
  pivot_longer(everything()) %>% mutate(value = round(value, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# ── plot: effects in H2 units, error band from the position's own MS_error ───
plot_df <- eff %>%
  pivot_longer(all_of(TERMS), names_to = "term", values_to = "H2") %>%
  mutate(chr = factor(chr, levels = chr_order), term = factor(term, levels = TERMS),
         pos_mb = pos / 1e6, se = sqrt(ms_error / 8)) %>%
  filter(!is.na(chr))

band <- plot_df %>% distinct(chr, pos_mb, se)

chr_lab_df <- plot_df %>% distinct(chr) %>% mutate(lab = chr_labels[as.character(chr)])

p <- ggplot(plot_df, aes(pos_mb, H2, colour = term)) +
  geom_ribbon(data = band, aes(x = pos_mb, ymin = -2*se, ymax = 2*se),
              fill = "grey80", alpha = 0.55, colour = NA, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ chr, ncol = 1, scales = "free_x") +
  geom_text(data = chr_lab_df, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.3, vjust = 1.4, size = 3, colour = "grey30",
            fontface = "bold", inherit.aes = FALSE) +
  scale_colour_manual(values = TERM_COLS, name = NULL) +
  scale_x_continuous(expand = expansion(0), minor_breaks = seq(0, 40, 1)) +
  labs(x = "Position (Mb)", y = expression("Cutler" ~ H^2 ~ "(deviation)"),
       subtitle = "grey band: +/- 2 SE from the within-cell odd-vs-even error at that position") +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 7, colour = "grey35"),
        strip.background = element_blank(), strip.text = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing = unit(0.4, "lines"))

dir.create("figures", showWarnings = FALSE)
png(FIG, width = 8, height = 6.5, units = "in", res = 150)
print(p)
invisible(dev.off())
cat("\nWrote:", FIG, "\n")
