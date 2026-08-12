# plot_H2_components.R — the heritability at each position, split into the
# component due to the overall mean, to sex, to diet, and to sex:diet.
#
# At each window the four cells are turned into the four orthonormal contrasts
# (mean, sex, diet, interaction), whose squares add to the total. But a squared
# contrast is noise-inflated -- E[x^2] = signal^2 + noise -- so at a background
# window all four squares are pure noise of similar size and the stack shows a
# tidy 25/25/25/25 that looks like a result and is not.
#
# Instead each component is the PRODUCT of the odd and even estimate:
#
#   component = contrast_odd * contrast_even        E[.] = signal^2
#
# Noise is independent between halves, so it cancels in expectation and
# background windows sit at zero rather than splitting into four equal parts.
# Individual products are noisy and can go negative, so they are averaged over a
# rolling window before plotting.
#
# Run from the XQTL2-dev repo ROOT, locally:
#   Rscript scripts_oneoffs/AGE_SY/splithalf/plot_H2_components.R [long_file.txt.gz]

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
IN   <- if (length(args) >= 1) args[1] else
  "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
FIG  <- if (length(args) >= 2) args[2] else "figures/AGE_SY_H2_components.png"

ROLL      <- 21L    # windows in the rolling mean (odd number; ~1.6 Mb at 75 kb)
TERMS     <- c("mean", "sex", "diet", "sex:diet")
TERM_COLS <- c(mean = "#7F7F7F", sex = "#1F78B4", diet = "#33A02C", `sex:diet` = "#E31A1C")
chr_order <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

if (!file.exists(IN)) stop("input not found: ", IN)
long <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble()

wide <- long %>%
  select(chr, pos, sex, sugar, half, Cutl_H2) %>%
  unite("cell", sugar, sex) %>%
  pivot_wider(names_from = cell, values_from = Cutl_H2) %>%
  drop_na(SY10_F, SY20_F, SY10_M, SY20_M) %>%
  transmute(chr, pos, half,
            mean       = (SY10_F + SY20_F + SY10_M + SY20_M) / 2,
            sex        = (SY10_F + SY20_F - SY10_M - SY20_M) / 2,
            diet       = (SY20_F + SY20_M - SY10_F - SY10_M) / 2,
            `sex:diet` = (SY10_F - SY20_F - SY10_M + SY20_M) / 2)

# ── odd x even, so the noise floor cancels ───────────────────────────────────
# Each contrast is centred on its own genome-wide level first. Without this the
# mean term is dominated by the fact that H2 is nonzero everywhere -- a constant
# baseline b makes the mean contrast a constant 2b whose product is a large
# positive number at every window, swamping the stack with the background level.
# Centring makes every component read as departure from the genome-wide norm, so
# background sits at zero and only real local structure shows.
comp <- wide %>%
  pivot_longer(all_of(TERMS), names_to = "term", values_to = "value") %>%
  pivot_wider(names_from = half, values_from = value) %>%
  drop_na(odd, even) %>%
  group_by(term) %>%
  mutate(odd = odd - mean(odd), even = even - mean(even)) %>%
  ungroup() %>%
  mutate(component = odd * even)

# Rolling mean with PARTIAL windows at the ends. Backfilling the edges with the
# raw unsmoothed values instead puts an unsmoothed point next to smoothed ones,
# which shows up as a spurious spike at every chromosome end.
roll_mean <- function(x, k) {
  n <- length(x)
  h <- (k - 1L) %/% 2L
  vapply(seq_len(n),
         function(i) mean(x[max(1L, i - h):min(n, i + h)], na.rm = TRUE),
         numeric(1))
}

comp <- comp %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  filter(!is.na(chr)) %>%
  arrange(chr, term, pos) %>%
  group_by(chr, term) %>%
  mutate(smooth = roll_mean(component, ROLL)) %>%
  ungroup() %>%
  mutate(pos_mb = pos / 1e6)

totals <- comp %>%
  group_by(chr, pos_mb) %>%
  summarise(total = sum(smooth), .groups = "drop")

# Shares are only meaningful where there is a total to take a share OF. A
# quantile cutoff is no good -- if signal covers little of the genome, a
# quantile just admits noise, and noise divided by noise fills the panel with
# vivid nonsense. Most windows are null, so the spread of the totals IS the
# noise scale: cut at median + 3 MAD, which excludes essentially all of it.
noise_scale <- mad(totals$total, na.rm = TRUE)
CUT <- max(median(totals$total, na.rm = TRUE) + 3 * noise_scale, 0)
kept <- mean(totals$total > CUT, na.rm = TRUE)
cat(sprintf("Shares drawn where smoothed total > %.3g (median + 3 MAD); %.1f%% of windows\n",
            CUT, 100 * kept))
if (kept < 0.002)
  cat("Very little of the genome clears the cut -- there may be no real structure.\n")

shares <- comp %>%
  left_join(totals, by = c("chr", "pos_mb")) %>%
  filter(total > CUT) %>%
  group_by(chr, pos_mb) %>%
  mutate(share = pmax(smooth, 0) / sum(pmax(smooth, 0))) %>%
  ungroup() %>%
  mutate(term = factor(term, levels = TERMS))

# geom_area interpolates straight across gaps, which would invent shares in the
# blank regions. Break each run of retained windows into its own group.
seg <- shares %>%
  distinct(chr, pos_mb) %>% arrange(chr, pos_mb) %>%
  group_by(chr) %>%
  mutate(seg = cumsum(c(TRUE, diff(pos_mb) > 1.5 * (ROLL * 0.075)))) %>%
  ungroup() %>%
  mutate(seg = paste0(chr, "_", seg))
shares <- shares %>% left_join(seg, by = c("chr", "pos_mb"))

p_total <- ggplot(totals, aes(pos_mb, total)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = CUT, colour = "grey60", linetype = "dashed", linewidth = 0.3) +
  geom_line(linewidth = 0.4) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  labs(y = "total signal", x = NULL) +
  theme_classic(base_size = 8) +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing = unit(0.15, "lines"))

p_share <- ggplot(shares, aes(pos_mb, share, fill = term, group = interaction(term, seg))) +
  # Force each facet to span the WHOLE chromosome, not just the retained
  # windows -- otherwise the kept regions stretch to fill the panel and the two
  # plots end up on different x scales, which reads as though signal is
  # everywhere.
  geom_blank(data = transform(totals, share = 0), aes(pos_mb, share),
             inherit.aes = FALSE) +
  geom_area(position = "stack", colour = NA) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = TERM_COLS, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  labs(x = "Position (Mb)", y = "share of signal") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_blank(), panel.spacing = unit(0.15, "lines"))

# Match the two plots' widths so the panels line up column for column; the
# y-axis labels differ in width and would otherwise shift the panels sideways.
g1 <- ggplot2::ggplotGrob(p_total)
g2 <- ggplot2::ggplotGrob(p_share)
maxw <- grid::unit.pmax(g1$widths, g2$widths)
g1$widths <- maxw
g2$widths <- maxw
g <- gridExtra::arrangeGrob(g1, g2, ncol = 1, heights = c(1, 1.6))
dir.create("figures", showWarnings = FALSE)
png(FIG, width = 11, height = 5.5, units = "in", res = 150)
grid::grid.draw(g)
invisible(dev.off())
cat("Wrote:", FIG, "\n")

cat("\nWhere each component peaks (smoothed odd x even):\n")
comp %>% group_by(term) %>% slice_max(smooth, n = 1) %>%
  transmute(term, chr, pos_mb = round(pos_mb, 2), signal = signif(smooth, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
