#!/usr/bin/env Rscript
###############################################################################
# anova_scan.R — the 2022 paper's test, run on the current pipeline's haplotypes.
#
# Macdonald & Long 2022 (G3, jkac279) tested, at each window,
#     ASF ~ H + TRT + H:TRT,   H:TRT against the R:H:TRT error term,
# with ASF = asin(sqrt(founder frequency)), H = 8 founders, R = 4 control/selected
# pairs (2 replicates x 2 sexes). Every replicate counts equally and the error is
# how much the replicates disagree. The pipeline's Wald test instead weights
# replicates by (effective) pool size and models the variance from sampling
# plus the lsei error, with no between-replicate term.
#
# Running the paper's statistic on the SAME smoothed frequencies the Wald test
# used separates "the test changed" from "the haplotype calls changed".
#
# Closed form, with d[r,h] = ASF_Z - ASF_C for replicate r and founder h:
#   F = [ R * sum_h (dbar_h - dbar)^2 / (H-1) ]
#     / [ sum_rh (d - dbar_h - dbar_r + dbar)^2 / ((R-1)(H-1)) ]
# Checked against anova(lm(y ~ R*H*TRT)) to 6 decimals.
#
# Also writes, for the two published peaks (Table 1: chr2R 12,294,089 and
# chr3L 6,162,469), the per-replicate founder shifts -- the pipeline's version
# of the paper's Fig. 4 -- and a one-line peak summary for every scan listed.
#
# Usage (from the repo root):
#   Rscript scripts_oneoffs/malathion/anova_scan.R \
#       --dir  process/malathion_2022 \
#       --scan malathion_2022_all \
#       --also malathion_2022_R1,malathion_2022_R2,malathion_2022_R3,malathion_2022_R4 \
#       --out  logs/malathion
###############################################################################
suppressPackageStartupMessages({ library(tidyverse); library(slider) })
options(width = 200)

args <- commandArgs(trailingOnly = TRUE)
p <- list(also = "")
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--dir"  = { p$dir  <- args[i+1]; i <- i+2L },
    "--scan" = { p$scan <- args[i+1]; i <- i+2L },
    "--also" = { p$also <- args[i+1]; i <- i+2L },
    "--out"  = { p$out  <- args[i+1]; i <- i+2L },
    stop("unknown argument: ", args[i]))
}
stopifnot(!is.null(p$dir), !is.null(p$scan), !is.null(p$out))
dir.create(p$out, showWarnings = FALSE, recursive = TRUE)

scan_dir <- function(s) file.path(p$dir, "Scans", s)
read_scan  <- function(s) as_tibble(read.table(file.path(scan_dir(s), paste0(s, ".scan.txt")), header = TRUE))
read_means <- function(s) as_tibble(read.table(file.path(scan_dir(s), paste0(s, ".meansBySample.txt")), header = TRUE))

# ── ANOVA at every window of the all-pools scan ──────────────────────────────
cat("reading", p$scan, "\n")
means <- read_means(p$scan)                       # chr pos TRT REP founder freq
wald  <- read_scan(p$scan) %>% select(chr, pos, cM, Wald_log10p)

d <- means %>%
  mutate(asf = asin(sqrt(pmin(pmax(freq, 0), 1)))) %>%
  select(chr, pos, TRT, REP, founder, asf) %>%
  pivot_wider(names_from = TRT, values_from = asf) %>%
  filter(!is.na(C), !is.na(Z)) %>%
  mutate(d = Z - C)

anova_window <- function(df) {
  m <- xtabs(d ~ REP + founder, df)               # R x H
  R <- nrow(m); H <- ncol(m)
  if (R < 2 || H < 2) return(tibble(F = NA_real_, anova_log10p = NA_real_, R = R, H = H))
  e  <- m - outer(rowMeans(m), colMeans(m), "+") + mean(m)
  Fv <- (R * sum((colMeans(m) - mean(m))^2) / (H - 1)) / (sum(e^2) / ((R - 1) * (H - 1)))
  tibble(F = Fv, anova_log10p = -pf(Fv, H - 1, (R - 1) * (H - 1), lower.tail = FALSE, log.p = TRUE) / log(10),
         R = R, H = H)
}

cat("ANOVA over", n_distinct(paste(d$chr, d$pos)), "windows\n")
an <- d %>% group_by(chr, pos) %>% group_modify(~ anova_window(.x)) %>% ungroup()
cat(sprintf("replicates per window: %s   founders: %s\n",
            paste(sort(unique(an$R)), collapse = "/"), paste(sort(unique(an$H)), collapse = "/")))

# concat_Chromosome_Scans.R slides the Wald over +/-5 windows; do the same here.
out <- wald %>% inner_join(an, by = c("chr", "pos")) %>%
  group_by(chr) %>% arrange(chr, pos) %>%
  mutate(anova_log10p_raw = anova_log10p,
         anova_log10p = slide_dbl(anova_log10p_raw, mean, .before = 5, .after = 5)) %>%
  ungroup()

f_tab <- file.path(p$out, paste0(p$scan, ".anova.txt.gz"))
write.table(out %>% select(chr, pos, cM, Wald_log10p, anova_log10p, anova_log10p_raw, F),
            gzfile(f_tab), row.names = FALSE, quote = FALSE)
cat("wrote", f_tab, "\n")

# ── Figure: Wald (pipeline) over ANOVA (paper), same haplotypes ──────────────
chr_lv <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
long <- out %>%
  select(chr, pos, Wald_log10p, anova_log10p) %>%
  pivot_longer(-c(chr, pos), names_to = "test", values_to = "lp") %>%
  mutate(chr = factor(chr, chr_lv),
         test = factor(test, c("Wald_log10p", "anova_log10p"),
                       c("Pipeline Wald test", "2022 paper ANOVA, same haplotypes")))
genes <- tibble(name = c("Cyp6g1", "Mdr65", "Ace"), chr = factor(c("chr2R", "chr3L", "chr3R"), chr_lv),
                pos = c(12185667, 6240500, 9069721))
g <- ggplot(long, aes(pos / 1e6, lp)) +
  geom_vline(data = genes, aes(xintercept = pos / 1e6), colour = "grey70", linetype = 2) +
  geom_line(linewidth = 0.3) +
  facet_grid(test ~ chr, scales = "free", space = "free_x", switch = "x") +
  geom_hline(yintercept = 4, colour = "red", linetype = 3) +
  labs(x = "Mb", y = expression(-log[10] * P), title = paste0(p$scan, ": the two tests on one set of haplotypes")) +
  theme_bw(base_size = 9) + theme(panel.spacing.x = unit(0, "pt"), strip.placement = "outside")
f_png <- file.path(p$out, paste0(p$scan, ".wald_vs_anova.png"))
ggsave(f_png, g, width = 10, height = 5, dpi = 200)
cat("wrote", f_png, "\n")

# ── Per-replicate founder shifts at the published peaks (the paper's Fig. 4) ──
peaks <- tibble(label = c("chr2R_Cyp6g1", "chr3L_Mdr65"), chr = c("chr2R", "chr3L"), peak = c(12294089, 6162469))
shift <- map_dfr(seq_len(nrow(peaks)), function(k) {
  w <- means %>% filter(chr == peaks$chr[k]) %>% distinct(pos) %>%
    mutate(dist = abs(pos - peaks$peak[k])) %>% slice_min(dist, n = 1) %>% pull(pos)
  means %>% filter(chr == peaks$chr[k], pos == w) %>%
    select(REP, TRT, founder, freq) %>%
    pivot_wider(names_from = TRT, values_from = freq) %>%
    mutate(peak = peaks$label[k], window = w, shift = Z - C) %>%
    select(peak, window, REP, founder, C, Z, shift)
})
shift_wide <- shift %>% select(peak, window, founder, REP, shift) %>%
  mutate(REP = paste0("R", REP)) %>%
  pivot_wider(names_from = REP, values_from = shift) %>%
  rowwise() %>% mutate(mean = mean(c_across(starts_with("R")))) %>% ungroup() %>%
  mutate(across(where(is.numeric) & !window, ~ round(.x, 3)))
f_shift <- file.path(p$out, paste0(p$scan, ".peak_shifts.txt"))
sink(f_shift)
cat("Founder frequency shift (selected - control) at the window nearest each published peak,\n")
cat("per replicate pair, from the smoothed haplotype frequencies the Wald test used.\n")
cat("REP labels as in the design file of the scan (design_2022: R1 = F rep1, R2 = M rep1, R3 = F rep2, R4 = M rep2).\n")
cat("Compare to Fig. 4 of Macdonald & Long 2022.\n\n")
print(as.data.frame(shift_wide), row.names = FALSE)
cat("\nSum of squared shifts per replicate (how much signal each pair carries):\n")
print(as.data.frame(shift %>% group_by(peak, REP) %>% summarise(sum_sq = round(sum(shift^2), 4), .groups = "drop") %>%
        pivot_wider(names_from = REP, values_from = sum_sq, names_prefix = "R")), row.names = FALSE)
sink()
cat("wrote", f_shift, "\n")

# ── Peak summary across every scan: max in the two intervals ─────────────────
scans <- c(p$scan, strsplit(p$also, ",")[[1]]) %>% discard(~ .x == "")
regions <- tibble(label = c("chr2R 10-14 Mb", "chr3L 5-7.5 Mb"), chr = c("chr2R", "chr3L"), lo = c(10e6, 5e6), hi = c(14e6, 7.5e6))
summ <- map_dfr(scans, function(s) {
  f <- file.path(scan_dir(s), paste0(s, ".scan.txt"))
  if (!file.exists(f)) return(tibble(scan = s, note = "scan.txt missing"))
  sc <- read_scan(s)
  map_dfr(seq_len(nrow(regions)), function(k) {
    r <- sc %>% filter(chr == regions$chr[k], pos >= regions$lo[k], pos <= regions$hi[k])
    tibble(scan = s, region = regions$label[k],
           max_wald = round(max(r$Wald_log10p, na.rm = TRUE), 2),
           at_Mb = round(r$pos[which.max(r$Wald_log10p)] / 1e6, 3))
  }) %>% pivot_wider(names_from = region, values_from = c(max_wald, at_Mb))
})
if ("anova_log10p" %in% names(out)) {
  summ <- bind_rows(summ, map_dfr(seq_len(nrow(regions)), function(k) {
    r <- out %>% filter(chr == regions$chr[k], pos >= regions$lo[k], pos <= regions$hi[k])
    tibble(scan = paste0(p$scan, " [paper ANOVA]"), region = regions$label[k],
           max_wald = round(max(r$anova_log10p, na.rm = TRUE), 2),
           at_Mb = round(r$pos[which.max(r$anova_log10p)] / 1e6, 3))
  }) %>% pivot_wider(names_from = region, values_from = c(max_wald, at_Mb)))
}
f_summ <- file.path(p$out, "peak_summary.txt")
sink(f_summ)
cat("Maximum -log10P in the two published QTL intervals, per scan.\n")
cat("Published (2022 ANOVA, 8 pools): chr2R ~5.2 at 12.29 Mb, chr3L ~4.9 at 6.16 Mb.\n")
cat("Stuart's rerun of the current pipeline: chr2R ~19.6, chr3L ~5.2.\n\n")
print(as.data.frame(summ), row.names = FALSE)
sink()
cat("wrote", f_summ, "\n")
cat(readLines(f_summ), sep = "\n")
