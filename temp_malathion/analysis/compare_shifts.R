#!/usr/bin/env Rscript
###############################################################################
# compare_shifts.R — within replicate 1, is the gap between the two pipelines
# in the haplotype CALLER or in the TEST?
#
# Same four pools, same flies, same reads, so nothing here is experimental
# sampling error. Either the two pipelines estimate different founder frequency
# shifts at the two peaks (the caller), or they estimate the same shifts and
# turn them into different p-values (the test). This prints both pipelines'
# shifts side by side at the published peaks.
#
# 2022 shifts are computed from the paper's published 1.5 cM frequencies
# (figshare 10.25387/g3.21200395); current-pipeline shifts are read from
# logs/malathion/malathion_rep1_all.peak_shifts.txt, written by the cluster run.
#
#   Rscript temp_malathion/analysis/compare_shifts.R
###############################################################################
suppressPackageStartupMessages({library(tidyverse); library(data.table)})

haplo_file <- "temp_malathion/data/malathion_haps_genetic_1p5cM.txt.gz"
shift_file <- "logs/malathion/malathion_rep1_all.peak_shifts.txt"
out_tab    <- "temp_malathion/analysis/shift_comparison.txt"
stopifnot(file.exists(haplo_file), file.exists(shift_file))

# published peaks (Table 1), in cM
peaks <- tibble(peak = c("chr2R", "chr3L"), CHROM = c("chr2R", "chr3L"), cM_peak = c(65.607, 16.627))

# ── 2022 pipeline: shifts for the rep-1 pools at the nearest window ─────────
hf <- fread(haplo_file, header = TRUE) %>% filter(pool %in% c("C1F","S1F","C1M","S1M"))
old <- map_dfr(seq_len(nrow(peaks)), function(k) {
  w <- hf %>% filter(chr == peaks$CHROM[k]) %>% distinct(cM) %>%
    slice_min(abs(cM - peaks$cM_peak[k]), n = 1) %>% pull(cM)
  hf %>% filter(chr == peaks$CHROM[k], cM == w) %>%
    mutate(TRT = str_sub(pool, 1, 1), sex = str_sub(pool, 3, 3)) %>%
    select(founder, sex, TRT, freq) %>%
    pivot_wider(names_from = TRT, values_from = freq) %>%
    transmute(peak = peaks$peak[k], pipeline = "2022", founder,
              REP = if_else(sex == "F", "R1", "R2"), shift = S - C)
})

# ── current pipeline: the table the cluster run wrote ───────────────────────
raw <- readLines(shift_file)
hdr <- grep("^\\s*peak\\s+window\\s+founder", raw)
new <- read.table(text = raw[hdr:(hdr + 16)], header = TRUE) %>%
  select(peak, founder, R1, R2) %>%
  mutate(peak = if_else(grepl("2R", peak), "chr2R", "chr3L")) %>%
  pivot_longer(c(R1, R2), names_to = "REP", values_to = "shift") %>%
  mutate(pipeline = "current")

both <- bind_rows(old, new) %>%
  pivot_wider(names_from = pipeline, values_from = shift) %>%
  arrange(peak, REP, founder)

sumsq <- both %>% group_by(peak, REP) %>%
  summarise(`2022` = sum(`2022`^2), current = sum(current^2), .groups = "drop")
ratio <- sumsq %>% group_by(REP) %>%
  summarise(`2022` = `2022`[peak == "chr2R"] / `2022`[peak == "chr3L"],
            current = current[peak == "chr2R"] / current[peak == "chr3L"], .groups = "drop")
agree <- both %>% group_by(peak) %>%
  summarise(r = cor(`2022`, current), .groups = "drop")

sink(out_tab)
cat("Founder frequency shift (selected - control) at the two published peaks,\n")
cat("replicate 1 only, the SAME four pools through each pipeline.\n")
cat("R1 = females, R2 = males.\n\n")
print(as.data.frame(both %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\nCorrelation between the two pipelines' shifts, across founders and sexes:\n")
print(as.data.frame(agree %>% mutate(r = round(r, 3))), row.names = FALSE)
cat("\nSum of squared shifts (the signal each test is given):\n")
print(as.data.frame(sumsq %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)
cat("\n2R / 3L ratio of that signal, per pipeline:\n")
print(as.data.frame(ratio %>% mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
cat("\nIf the ratios agree, both callers see the same relative signal at the two\n")
cat("peaks and the difference in peak HEIGHT is the test. If they disagree, the\n")
cat("caller is contributing too.\n")
sink()
cat(readLines(out_tab), sep = "\n")

# ── what each test does with those identical shifts ─────────────────────────
# Appended: the paper's ANOVA table at each peak, replicate 1 only, so the
# numerator (the treatment-by-founder signal) and the denominator (the yardstick,
# which absorbs how much the two sexes disagree) can be read separately.
cat("\n\n")
tab2 <- map_dfr(seq_len(nrow(peaks)), function(k) {
  w <- hf %>% filter(chr == peaks$CHROM[k]) %>% distinct(cM) %>%
    slice_min(abs(cM - peaks$cM_peak[k]), n = 1) %>% pull(cM)
  d <- hf %>% filter(chr == peaks$CHROM[k], cM == w) %>%
    mutate(TRT = as.factor(str_sub(pool, 1, 1)),
           HAP = as.factor(founder), asf = asin(sqrt(freq)))
  tt <- tapply(d$asf, d$HAP, mean); keep <- names(tt)[tt > 0.14]
  d3 <- d %>% filter(HAP %in% keep) %>% droplevels()
  a  <- anova(lm(asf ~ TRT + HAP + TRT:HAP, data = d3))
  tibble(peak = peaks$peak[k], founders_kept = nlevels(d3$HAP),
         MS_signal = a[3, 3], df_signal = a[3, 1],
         MS_error = a[4, 3], df_error = a[4, 1],
         F = a[3, 3] / a[4, 3],
         mlog10p = -pf(a[3,3]/a[4,3], a[3,1], a[4,1], lower.tail = FALSE, log.p = TRUE)/log(10))
})
sink(out_tab, append = TRUE)
cat("\n\nThe paper's ANOVA at each peak, replicate 1 only, from the published frequencies.\n")
cat("MS_signal is the treatment-by-founder mean square; MS_error is the residual,\n")
cat("which is where disagreement between the two sexes goes.\n\n")
print(as.data.frame(tab2 %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
sink()
cat(readLines(out_tab), sep = "\n")
