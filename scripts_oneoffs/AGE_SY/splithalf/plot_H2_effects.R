# plot_H2_effects.R — the same 5-panel layout as the 4-scan figure, but the
# four lines are the COMPONENTS of Cutler H2 rather than the four treatments.
#
# At each window the four treatment values are rewritten as
#
#   H2(sugar, sex) = mean  +/- sex  +/- sugar  +/- interaction
#
#   mean         average of the four                       (a level)
#   sex          (females - males) / 2                     (a deviation)
#   sugar        (SY20 - SY10) / 2                         (a deviation)
#   interaction  (SY10_F - SY20_F - SY10_M + SY20_M) / 4   (a deviation)
#
# so at the chr3L QTL, where the four treatments read about 2.8, 2.8, 4.9, 5.5,
# the lines read mean 4.0, sex -1.2, sugar 0.15, interaction 0.15. Everything is
# in H2 units and adds back to the four treatment values exactly.
#
# The two halves are averaged for display. Set SHOW_HALVES to draw them
# separately instead -- where the two traces agree, the feature is real.
#
# Run from the XQTL2-dev repo ROOT, locally:
#   Rscript scripts_oneoffs/AGE_SY/splithalf/plot_H2_effects.R [long_file.txt.gz]

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
IN   <- if (length(args) >= 1) args[1] else
  "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
FIG  <- if (length(args) >= 2) args[2] else "figures/AGE_SY_H2_effects.png"
SHOW_HALVES <- FALSE

TERMS     <- c("mean", "sex", "sugar", "interaction")
TERM_COLS <- c(mean = "black", sex = "#1F78B4", sugar = "#33A02C",
               interaction = "#E31A1C")
chr_order  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
chr_labels <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")

# dm6 euchromatin boundaries (Huynh et al. 2023 PLoS Genet 19:e1010439, Table S2)
HET_BOUNDS <- tribble(
  ~chr,    ~eu_start, ~eu_end,
  "chrX",   2.5,      21.2,
  "chr2L",  0.5,      22.9,
  "chr2R",  1.3,      25.1,
  "chr3L",  0.7,      24.0,
  "chr3R",  4.5,      32.0
)

if (!file.exists(IN)) stop("input not found: ", IN)
long <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble()

cells <- long %>%
  select(chr, pos, sex, sugar, half, Cutl_H2) %>%
  unite("cell", sugar, sex) %>%
  pivot_wider(names_from = cell, values_from = Cutl_H2) %>%
  drop_na(SY10_F, SY20_F, SY10_M, SY20_M)

if (!SHOW_HALVES) {
  cells <- cells %>%
    group_by(chr, pos) %>%
    summarise(across(c(SY10_F, SY20_F, SY10_M, SY20_M), mean), .groups = "drop") %>%
    mutate(half = "both")
}

eff <- cells %>%
  transmute(chr, pos, half,
            mean        = (SY10_F + SY20_F + SY10_M + SY20_M) / 4,
            sex         = (SY10_F + SY20_F - SY10_M - SY20_M) / 4 * 2,
            sugar       = (SY20_F + SY20_M - SY10_F - SY10_M) / 4 * 2,
            interaction = (SY10_F - SY20_F - SY10_M + SY20_M) / 4) %>%
  pivot_longer(all_of(TERMS), names_to = "term", values_to = "H2") %>%
  mutate(chr = factor(chr, levels = chr_order),
         term = factor(term, levels = TERMS),
         pos_mb = pos / 1e6) %>%
  filter(!is.na(chr))

het_rects <- HET_BOUNDS %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  left_join(eff %>% group_by(chr) %>% summarise(xmax = max(pos_mb), .groups = "drop"),
            by = "chr") %>%
  { bind_rows(transmute(., chr, xmin = 0,      xmax2 = eu_start),
              transmute(., chr, xmin = eu_end, xmax2 = xmax)) } %>%
  rename(xmax = xmax2)

chr_lab_df <- eff %>% distinct(chr) %>% mutate(lab = chr_labels[as.character(chr)])

p <- ggplot(eff, aes(pos_mb, H2, colour = term)) +
  geom_rect(data = het_rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.5, colour = NA, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  { if (SHOW_HALVES) geom_line(aes(linetype = half), linewidth = 0.5)
    else geom_line(linewidth = 0.5) } +
  facet_wrap(~ chr, ncol = 1, scales = "free_x") +
  geom_text(data = chr_lab_df, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.3, vjust = 1.4, size = 3, colour = "grey30",
            fontface = "bold", inherit.aes = FALSE) +
  scale_colour_manual(values = TERM_COLS, name = NULL) +
  scale_x_continuous(expand = expansion(0), minor_breaks = seq(0, 40, 1)) +
  labs(x = "Position (Mb)", y = expression("Cutler" ~ H^2)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top",
        strip.background = element_blank(), strip.text = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing = unit(0.4, "lines"))

dir.create("figures", showWarnings = FALSE)
png(FIG, width = 8, height = 6.5, units = "in", res = 150)
print(p)
invisible(dev.off())
cat("Wrote:", FIG, "\n")

cat("\nGenome-wide median of each component (H2 units):\n")
eff %>% group_by(term) %>% summarise(median = round(median(H2), 3), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nLargest absolute deviation per component, and where:\n")
eff %>% group_by(term) %>% slice_max(abs(H2), n = 1) %>%
  transmute(term, chr, pos_mb = round(pos_mb, 2), H2 = round(H2, 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)
