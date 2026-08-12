# anova_H2.R — the nested partition of Cutler H2 at each position.
#
#   log2 H2  ~  sex + diet + sex:diet + rep %in% (sex:diet)
#
# H2 IS A VARIANCE, so the four treatments are compared as RATIOS, not
# differences: the model is additive on the log scale and the effects are fold
# changes in heritability. Two things in the data say this is the right scale:
#
#   * no value is <= 0 (minimum 0.099), so nothing is dropped by taking logs
#   * the within-cell spread scales with the level -- SD runs 0.130 to 0.449 as
#     H2 runs 0.46 to 1.64, while the coefficient of variation stays flat at
#     0.28-0.32. Constant CV means constant variance on the log scale, i.e. the
#     error is multiplicative
#
# The odd and even halves are the replicates: two independent estimates per cell
# per window, so the four within-cell odd-vs-even differences give a 4 df error
# term at every position, and no genome-wide noise model is needed.
#
# With the four cell values on log2 as a = SY10_F, b = SY20_F, c = SY10_M,
# d = SY20_M:
#
#   sex       (a+b)/2 - (c+d)/2     log2 fold change, females / males
#   diet      (b+d)/2 - (a+c)/2     log2 fold change, SY20 / SY10
#   sex:diet  (b-a) - (d-c)         how much the diet fold change differs
#                                   between the sexes
#
# 0 means no effect, -1 means twofold lower, +1 twofold higher.
#
# Each cell mean averages n=2 halves, so var(cell mean) = MS_error / 2 and each
# contrast carries its own SE from its coefficients: sex and diet sqrt(MS/2),
# the interaction sqrt(2*MS).
#
# MS_error on 4 df is far too noisy per window, and windows step 5 kb under a
# 75 kb haplotype window so neighbours are near-duplicates. It is pooled over
# POOL_BP. That pools the ERROR, not the effects.
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript scripts_oneoffs/AGE_SY/splithalf/anova_H2.R [long_file.txt.gz]

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
IN   <- if (length(args) >= 1) args[1] else
  "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
OUT  <- "process/AGE_SY_splithalf/H2_anova_by_window.txt.gz"
FIG  <- "figures/AGE_SY_H2_anova.png"
POOL_BP <- 5e5

TERMS     <- c("sex", "diet", "sex:diet")
TERM_COLS <- c(sex = "#1F78B4", diet = "#33A02C", `sex:diet` = "#E31A1C")
chr_order  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
chr_labels <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")

if (!file.exists(IN)) stop("input not found: ", IN)
long <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble()

if (any(long$Cutl_H2 <= 0, na.rm = TRUE))
  stop(sum(long$Cutl_H2 <= 0, na.rm = TRUE),
       " values are <= 0; the log scale is not usable without handling them")

w <- long %>%
  mutate(l = log2(Cutl_H2)) %>%
  select(chr, pos, sex, sugar, half, l) %>%
  unite("cell", sugar, sex) %>%
  pivot_wider(names_from = c(cell, half), values_from = l, names_glue = "{cell}_{half}") %>%
  drop_na()

cm <- w %>% transmute(
  chr, pos,
  a = (SY10_F_odd + SY10_F_even) / 2,     # SY10 female
  b = (SY20_F_odd + SY20_F_even) / 2,     # SY20 female
  c = (SY10_M_odd + SY10_M_even) / 2,     # SY10 male
  d = (SY20_M_odd + SY20_M_even) / 2,     # SY20 male
  ms_error = ((SY10_F_odd - SY10_F_even)^2 + (SY20_F_odd - SY20_F_even)^2 +
              (SY10_M_odd - SY10_M_even)^2 + (SY20_M_odd - SY20_M_even)^2) / 2 / 4)

cm <- cm %>%
  mutate(region = paste0(chr, "_", pos %/% POOL_BP)) %>%
  group_by(region) %>% mutate(ms_error = mean(ms_error), n_pool = n()) %>% ungroup()

eff <- cm %>%
  transmute(chr, pos, ms_error,
            level_H2   = 2^((a + b + c + d) / 4),        # geometric mean H2
            sex        = (a + b) / 2 - (c + d) / 2,
            diet       = (b + d) / 2 - (a + c) / 2,
            `sex:diet` = (b - a) - (d - c),
            se_sex     = sqrt(ms_error / 2),
            se_diet    = sqrt(ms_error / 2),
            `se_sex:diet` = sqrt(2 * ms_error))

write.table(eff, gzfile(OUT), row.names = FALSE, quote = FALSE, sep = "\t")
cat("Wrote:", OUT, "\n")
cat(sprintf("MS_error pooled over %.0f kb (%.0f windows, 4 df each)\n",
            POOL_BP / 1e3, median(cm$n_pool)))

fold <- function(x) ifelse(x >= 0, sprintf("%.2fx up", 2^x), sprintf("%.2fx down", 2^(-x)))

cat("\nGenome-wide medians (log2 fold change, and as a ratio):\n")
tibble(term = TERMS,
       log2 = c(median(eff$sex), median(eff$diet), median(eff$`sex:diet`))) %>%
  mutate(ratio = fold(log2), log2 = round(log2, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nAt the chr3L peak (9.3-9.5 Mb):\n")
pk <- eff %>% filter(chr == "chr3L", pos > 9.3e6, pos < 9.5e6)
tibble(term = TERMS,
       log2 = c(mean(pk$sex), mean(pk$diet), mean(pk$`sex:diet`)),
       se   = c(mean(pk$se_sex), mean(pk$se_diet), mean(pk$`se_sex:diet`))) %>%
  mutate(ratio = fold(log2), `SE (log2)` = round(se, 3),
         `n SE` = round(abs(log2 / se), 1), log2 = round(log2, 3)) %>%
  select(-se) %>% as.data.frame() %>% print(row.names = FALSE)
cat(sprintf("  geometric mean H2 there: %.2f\n", mean(pk$level_H2)))

plot_df <- eff %>%
  pivot_longer(all_of(TERMS), names_to = "term", values_to = "log2fc") %>%
  mutate(se = case_when(term == "sex" ~ se_sex, term == "diet" ~ se_diet,
                        TRUE ~ `se_sex:diet`),
         chr = factor(chr, levels = chr_order), term = factor(term, levels = TERMS),
         pos_mb = pos / 1e6) %>%
  filter(!is.na(chr))

# One band per position: use the sex/diet SE, and note the interaction's is 2x.
band <- plot_df %>% filter(term == "sex") %>% distinct(chr, pos_mb, se)
chr_lab_df <- plot_df %>% distinct(chr) %>% mutate(lab = chr_labels[as.character(chr)])

p <- ggplot(plot_df, aes(pos_mb, log2fc, colour = term)) +
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
  scale_y_continuous(sec.axis = sec_axis(~ 2^., name = "fold change",
                                         breaks = c(0.25, 0.5, 1, 2, 4))) +
  labs(x = "Position (Mb)", y = expression(log[2] ~ "fold change in Cutler" ~ H^2),
       title = "Is the heritability at this position different between sexes, between diets, or both",
       subtitle = paste("H2 is a variance, so treatments are compared as ratios. 0 = no difference.",
                        "Grey band is +/- 2 SE from the odd-vs-even difference within each treatment;",
                        "the sex:diet SE is 2x the band.")) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top",
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        strip.background = element_blank(), strip.text = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing = unit(0.4, "lines"))

dir.create("figures", showWarnings = FALSE)
png(FIG, width = 8.5, height = 6.5, units = "in", res = 150)
print(p)
invisible(dev.off())
cat("\nWrote:", FIG, "\n")
