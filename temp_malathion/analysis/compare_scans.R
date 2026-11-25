#!/usr/bin/env Rscript
###############################################################################
# compare_scans.R — why the malathion scan changed: separate the three things
# that differ between the 2022 figure and the current pipeline's output.
#
# The published scan and Stuart's rerun differ in THREE ways at once: the set of
# pools (8 vs 4), the haplotype calls (2022 lsei at 1.5 cM vs the current
# pipeline), and the test (the paper's ANOVA vs the pipeline's Wald). Changing
# all three at once explains nothing, so each panel changes one thing.
#
#   A  paper test, paper haplotypes, 8 pools    = the published Fig. 2
#   B  paper test, paper haplotypes, rep 1 only   A->B isolates dropping rep 2
#   C  Wald test,  new haplotypes,  rep 1 only    B->C is new haps + new test
#
# A and B are computed here with the paper's own ANOVA_test_function, taken
# verbatim from the figshare deposit (malathion_analysis.zip,
# get_haplotype_frequencies.txt), run on the paper's own published haplotype
# frequencies (malathion_haps_genetic_1p5cM.txt.gz, same deposit). So A is a
# reproduction of the published figure, not a redrawing of it.
#
# C comes from logs/malathion/malathion_rep1_all.anova.txt.gz, written by the
# current pipeline on the four rep-1 pools HPC3 holds.
#
# Usage, from the repo root:
#   Rscript temp_malathion/analysis/compare_scans.R <dir with the figshare file>
###############################################################################
suppressPackageStartupMessages({library(tidyverse); library(data.table); library(patchwork)})

args <- commandArgs(trailingOnly = TRUE)
fs_dir <- if (length(args)) args[1] else "temp_malathion/data"
haplo_file <- file.path(fs_dir, "malathion_haps_genetic_1p5cM.txt.gz")
new_file   <- "logs/malathion/malathion_rep1_all.anova.txt.gz"
out_png    <- "temp_malathion/analysis/scan_comparison.png"
out_tab    <- "temp_malathion/analysis/scan_comparison.txt"

# The paper's haplotype frequencies and its analysis code are the figshare
# supplement, doi 10.25387/g3.21200395. Fetched once, not committed.
if (!file.exists(haplo_file)) {
  dir.create(fs_dir, showWarnings = FALSE, recursive = TRUE)
  cat("fetching the published haplotype frequencies from figshare...\n")
  download.file("https://ndownloader.figshare.com/files/37586293", haplo_file, quiet = TRUE, mode = "wb")
}
stopifnot(file.exists(haplo_file), file.exists(new_file))

# ── the paper's test, verbatim from the figshare deposit ────────────────────
# (malathion_analysis/get_haplotype_frequencies.txt, lines 209-231)
ANOVA_test_function = function(df) {
	df = data.frame(df)
	# drop rare haplotypes, roughly below 2%: asin(sqrt(0.02)) = 0.1418971
	tt = tapply(df$asf,df$HAP,mean)
	tt2 = names(tt)[tt > 0.14]
	df3 = df %>% filter(HAP %in% tt2) %>% droplevels()
	out=anova(lm(asf~TRT + HAP + TRT:HAP, data=df3))
	F=out[3,3]/out[4,3]
	df1=out[3,1]
	df2=out[4,1]
	-pf(F,df1,df2,lower.tail=FALSE,log.p=TRUE)/log(10)
}

# ── the paper's input preparation, same source ───────────────────────────────
cat("reading", haplo_file, "\n")
haplo_freqs <- fread(haplo_file, header = TRUE)
trans <- haplo_freqs %>%
  mutate(TRT = str_sub(pool, 1, 1),
         REP = as.numeric(as.factor(str_sub(pool, 2, 3))),
         HAP = gsub("B", "", gsub("A", "HAP", founder)),
         rep_expt = str_sub(pool, 2, 2),
         HAP = as.factor(HAP), TRT = as.factor(TRT),
         asf = asin(sqrt(freq))) %>%
  drop_na() %>%
  select(CHROM = chr, POS = pos, cM, TRT, REP, HAP, asf, rep_expt)

run_scan <- function(df, label) {
  cat("scanning", label, ":", n_distinct(paste(df$CHROM, df$POS)), "windows,",
      n_distinct(df$REP), "pairs\n")
  df %>% mutate(REP = as.factor(REP)) %>%
    group_by(CHROM, POS, cM) %>% nest() %>%
    mutate(mlog10p = map_dbl(data, ~ ANOVA_test_function(.))) %>%
    select(-data) %>% ungroup() %>% mutate(panel = label)
}

A <- run_scan(trans, "A")                                  # all 8 pools
B <- run_scan(trans %>% filter(rep_expt == "1"), "B")       # rep 1 only

# ── the current pipeline on the four rep-1 pools ────────────────────────────
C <- read.table(new_file, header = TRUE) %>%
  transmute(CHROM = chr, POS = pos, cM, mlog10p = Wald_log10p, panel = "C")

# ── cumulative cM axis, chromosomes in the paper's order ────────────────────
chr_lv <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
all3 <- bind_rows(A, B, C) %>% filter(CHROM %in% chr_lv) %>%
  mutate(CHROM = factor(CHROM, chr_lv))
chrlen <- all3 %>% group_by(CHROM) %>% summarise(len = max(cM), .groups = "drop") %>%
  arrange(CHROM) %>% mutate(off = cumsum(lag(len, default = 0)) + 2 * (as.integer(CHROM) - 1))
all3 <- all3 %>% left_join(chrlen, by = "CHROM") %>% mutate(cMcum = cM + off)
mid  <- all3 %>% group_by(CHROM) %>% summarise(m = mean(range(cMcum)), .groups = "drop")

genes <- tibble(name = c("Cyp6g1", "Mdr65", "Ace"), CHROM = factor(c("chr2R","chr3L","chr3R"), chr_lv),
                cM = c(65.607, 16.627, 30.6)) %>%
  left_join(chrlen, by = "CHROM") %>% mutate(cMcum = cM + off)

lab <- c(A = "A.  2022 test, 2022 haplotypes, all 8 pools  (the published Fig. 2)",
         B = "B.  2022 test, 2022 haplotypes, replicate 1 only (4 pools)",
         C = "C.  new pipeline Wald test, new haplotypes, replicate 1 only (4 pools)")

p <- ggplot(all3, aes(cMcum, mlog10p)) +
  geom_vline(data = genes, aes(xintercept = cMcum), colour = "grey75", linetype = 2, linewidth = .3) +
  geom_hline(yintercept = 4, colour = "firebrick", linetype = 3, linewidth = .3) +
  geom_line(aes(colour = CHROM), linewidth = .35) +
  geom_text(data = genes, aes(x = cMcum, y = Inf, label = name), inherit.aes = FALSE,
            vjust = 1.4, size = 2.5, colour = "grey30") +
  facet_wrap(~ panel, ncol = 1, scales = "free_y",
             labeller = labeller(panel = lab)) +
  scale_x_continuous(breaks = mid$m, labels = sub("^chr", "", levels(mid$CHROM))) +
  scale_colour_manual(values = rep(c("#1f4e79", "#2e8b57"), 3)[1:5], guide = "none") +
  labs(x = "Position (cM)", y = expression(-log[10] * P),
       title = "Malathion X-QTL: what changed between the 2022 scan and the current pipeline",
       subtitle = "One difference at a time. A to B drops replicate 2; B to C changes both the haplotype caller and the test.") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(hjust = 0, face = "bold", size = 8),
        strip.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10))

ggsave(out_png, p, width = 9, height = 7, dpi = 200)
cat("wrote", out_png, "\n")

# ── the two published peaks, every panel ────────────────────────────────────
peaks <- tibble(peak = c("chr2R (Cyp6g1)", "chr3L (Mdr65)"),
                CHROM = factor(c("chr2R", "chr3L"), chr_lv), lo = c(55, 10), hi = c(75, 25))
tab <- map_dfr(seq_len(nrow(peaks)), function(k) {
  all3 %>% filter(CHROM == peaks$CHROM[k], cM >= peaks$lo[k], cM <= peaks$hi[k]) %>%
    group_by(panel) %>%
    summarise(peak = peaks$peak[k], max = round(max(mlog10p, na.rm = TRUE), 2),
              at_cM = round(cM[which.max(mlog10p)], 2), .groups = "drop")
}) %>% select(panel, peak, max, at_cM) %>%
  pivot_wider(names_from = peak, values_from = c(max, at_cM)) %>% arrange(panel)

sink(out_tab)
cat("Maximum -log10P at the two published QTL, one difference at a time.\n\n")
for (k in names(lab)) cat(lab[[k]], "\n")
cat("D.  new pipeline Wald test, new haplotypes, all 8 pools -- Stuart's rerun, approx 19.6 and 5.2\n\n")
print(as.data.frame(tab), row.names = FALSE)
cat("\nPublished Table 1 peaks: chr2R 12,294,089 bp (65.6 cM), chr3L 6,162,469 bp (16.6 cM).\n")
cat("A is computed here from the paper's own code and data, so it is a reproduction\n")
cat("of the published figure rather than a reading off it.\n")
sink()
cat(readLines(out_tab), sep = "\n")
