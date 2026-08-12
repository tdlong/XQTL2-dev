# partition_H2.R — split the spatial variation in Cutler H2 into main effects
# and interaction, with the noise measured rather than assumed.
#
# Reads the long file from gather_splithalf_H2.R and, at each window, forms the
# four orthonormal contrasts of the 2x2 (sex x sugar), separately in the odd and
# the even half:
#
#   mean   ( SY10_F + SY20_F + SY10_M + SY20_M ) / 2
#   sex    ( SY10_F + SY20_F - SY10_M - SY20_M ) / 2
#   sugar  ( SY20_F + SY20_M - SY10_F - SY10_M ) / 2
#   inter  ( SY10_F - SY20_F - SY10_M + SY20_M ) / 2
#
# Coefficients are +-1/2, so the four contrasts are orthonormal and their
# squares sum to the squares of the four cells: a real partition.
#
# The noise: it is independent between odd and even, so across windows
#   cov(contrast_odd, contrast_even)  estimates the SIGNAL variance
#   mean(var_odd, var_even)           is signal + noise
# The gap between them is noise. This needs no assumption that the four
# contrasts carry equal noise -- they do not, since SY10 and SY20 share their
# controls, so control noise cancels in sugar and interaction but not in sex.
#
# Variance is taken ACROSS WINDOWS, so a cell that is uniformly higher drops out
# and what is partitioned is spatial variation, not level. Level differences are
# reported separately.
#
# Run from the XQTL2-dev repo ROOT, locally:
#   Rscript scripts_oneoffs/AGE_SY/splithalf/partition_H2.R
#
# Optional first argument: path to the long file.

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
IN   <- if (length(args) >= 1) args[1] else
  "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
OUTDIR <- "process/AGE_SY_splithalf"
FIG    <- "figures/AGE_SY_H2_partition.png"

if (!file.exists(IN)) stop("input not found: ", IN)
long <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble()

# ── one row per window per half, four cells wide ─────────────────────────────
wide <- long %>%
  select(chr, pos, sex, sugar, half, Cutl_H2) %>%
  unite("cell", sugar, sex) %>%              # SY10_F, SY20_F, SY10_M, SY20_M
  pivot_wider(names_from = cell, values_from = Cutl_H2)

need <- c("SY10_F", "SY20_F", "SY10_M", "SY20_M")
if (!all(need %in% names(wide)))
  stop("missing cells: ", paste(setdiff(need, names(wide)), collapse = ", "))

wide <- wide %>% drop_na(all_of(need))

contrasts_df <- wide %>%
  transmute(chr, pos, half,
            mean  = (SY10_F + SY20_F + SY10_M + SY20_M) / 2,
            sex   = (SY10_F + SY20_F - SY10_M - SY20_M) / 2,
            sugar = (SY20_F + SY20_M - SY10_F - SY10_M) / 2,
            inter = (SY10_F - SY20_F - SY10_M + SY20_M) / 2)

TERMS <- c("mean", "sex", "sugar", "inter")

# ── odd vs even, window by window ────────────────────────────────────────────
paired <- contrasts_df %>%
  pivot_longer(all_of(TERMS), names_to = "term", values_to = "value") %>%
  pivot_wider(names_from = half, values_from = value) %>%
  drop_na(odd, even)

partition <- paired %>%
  group_by(term) %>%
  summarise(
    n_windows = n(),
    var_total = (var(odd) + var(even)) / 2,   # signal + noise
    var_signal = cov(odd, even),              # noise is independent -> signal
    reliability = cov(odd, even) / sqrt(var(odd) * var(even)),
    .groups = "drop"
  ) %>%
  mutate(
    var_noise = var_total - var_signal,
    pct_of_signal = 100 * pmax(var_signal, 0) / sum(pmax(var_signal, 0)),
    pct_of_total  = 100 * var_total / sum(var_total),
    term = factor(term, levels = TERMS)
  ) %>%
  arrange(term)

# ── block bootstrap ──────────────────────────────────────────────────────────
# Adjacent windows are not independent -- haplotypes are smoothed at 100 kb and
# linkage correlates neighbours -- so naive standard errors would be far too
# small. Resample contiguous blocks of windows within chromosome instead.
# Block length is set in BASE PAIRS, not in windows: windows slide in 5 kb steps
# while the haplotype window is 75 kb and the smoothing 100 kb, so neighbours
# overlap heavily and a block counted in windows would be far shorter than the
# correlation length -- giving intervals that are much too narrow.
BLOCK_BP <- 2e6
B        <- 400L

step_bp <- paired %>% distinct(chr, pos) %>% arrange(chr, pos) %>%
  group_by(chr) %>% summarise(s = median(diff(pos)), .groups = "drop") %>%
  pull(s) %>% median()
cat(sprintf("Window step %.0f kb; bootstrap blocks of %.1f Mb (%.0f windows)\n",
            step_bp / 1e3, BLOCK_BP / 1e6, BLOCK_BP / step_bp))

boot_ids <- paired %>%
  distinct(chr, pos) %>% arrange(chr, pos) %>%
  mutate(block = paste0(chr, "_", pos %/% BLOCK_BP))

paired_b <- paired %>% left_join(boot_ids, by = c("chr", "pos"))
blocks   <- unique(paired_b$block)

set.seed(1)
boot <- map_dfr(seq_len(B), function(b) {
  drawn <- sample(blocks, length(blocks), replace = TRUE)
  idx   <- unlist(lapply(drawn, function(bl) which(paired_b$block == bl)), use.names = FALSE)
  paired_b[idx, ] %>%
    group_by(term) %>%
    summarise(var_signal = cov(odd, even),
              reliability = cov(odd, even) / sqrt(var(odd) * var(even)),
              .groups = "drop") %>%
    mutate(pct_of_signal = 100 * pmax(var_signal, 0) / sum(pmax(var_signal, 0)),
           rep = b)
})

ci <- boot %>%
  group_by(term) %>%
  summarise(rel_lo = quantile(reliability, 0.025, na.rm = TRUE),
            rel_hi = quantile(reliability, 0.975, na.rm = TRUE),
            sig_lo = quantile(var_signal, 0.025, na.rm = TRUE),
            pct_lo = quantile(pct_of_signal, 0.025, na.rm = TRUE),
            pct_hi = quantile(pct_of_signal, 0.975, na.rm = TRUE),
            .groups = "drop")

partition <- partition %>% left_join(ci, by = "term") %>%
  mutate(real = sig_lo > 0)   # signal variance clearly above zero

cat("\nSpatial variation in Cutler H2, partitioned\n")
cat("(variance across windows; a uniformly higher cell contributes nothing)\n\n")
partition %>%
  transmute(term,
            `var signal` = signif(var_signal, 3),
            `var noise`  = signif(var_noise, 3),
            `reliability` = sprintf("%.3f [%.3f, %.3f]", reliability, rel_lo, rel_hi),
            `% of signal` = sprintf("%.1f [%.0f, %.0f]", pct_of_signal, pct_lo, pct_hi),
            `signal > 0` = ifelse(real, "yes", "no"),
            `% of raw var` = round(pct_of_total, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nreliability = cor(odd, even) for that contrast: the share of its spread\n")
cat("that reproduces across independent halves. Near 0 means noise.\n")
cat(sprintf("Intervals are 95%% from a moving-block bootstrap (%d blocks of %.1f Mb\n",
            length(blocks), BLOCK_BP/1e6))
cat("windows, resampled within chromosome), which respects the fact that\n")
cat("neighbouring windows are correlated.\n")
if (any(!partition$real))
  cat(sprintf("signal > 0 is 'no' for: %s. Read those shares as shares of\nsomething not distinguishable from nothing.\n",
              paste(partition$term[!partition$real], collapse = ", ")))
cat("A wide interval on a term whose signal IS real usually means it rests on\n")
cat("one or two loci -- resamples that miss them find nothing. Worth knowing.\n")

# ── level differences, which the variance partition deliberately ignores ─────
cat("\nGenome-wide level of each contrast (NOT part of the partition above):\n")
contrasts_df %>%
  pivot_longer(all_of(TERMS), names_to = "term", values_to = "value") %>%
  group_by(term, half) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = half, values_from = mean) %>%
  mutate(term = factor(term, levels = TERMS)) %>% arrange(term) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("A sugar level difference is confounded with SY20 selecting a larger\n")
cat("fraction than SY10 -- systematic, so it sits in both halves.\n")

write.table(partition, file.path(OUTDIR, "H2_partition_summary.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

# ── the components along the genome, odd and even overlaid ───────────────────
chr_order <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
plot_df <- contrasts_df %>%
  pivot_longer(all_of(TERMS), names_to = "term", values_to = "value") %>%
  mutate(chr = factor(chr, levels = chr_order),
         term = factor(term, levels = TERMS),
         pos_mb = pos / 1e6) %>%
  filter(!is.na(chr))

p <- ggplot(plot_df, aes(pos_mb, value, colour = half)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.35, alpha = 0.85) +
  facet_grid(term ~ chr, scales = "free", space = "free_x", switch = "y") +
  scale_colour_manual(values = c(odd = "#1F78B4", even = "#E31A1C"), name = NULL) +
  labs(x = "Position (Mb)", y = expression("Cutler" ~ H^2 ~ "contrast")) +
  theme_classic(base_size = 8) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.25, "lines"))

dir.create("figures", showWarnings = FALSE)
png(FIG, width = 11, height = 6.5, units = "in", res = 150)
print(p)
invisible(dev.off())
cat("\nWrote:", file.path(OUTDIR, "H2_partition_summary.txt"), "\n")
cat("Wrote:", FIG, "\n")
cat("\nWhere the two halves trace the same shape, that shape is real.\n")
