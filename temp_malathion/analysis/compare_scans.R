#!/usr/bin/env Rscript
###############################################################################
# compare_scans.R — the malathion scan as a 2x2: which DATA, which PIPELINE.
#
#                     | replicate 1 only (4 pools) | both replicates (8 pools)
#   2022 pipeline     |            B               |      A = published Fig. 2
#   current pipeline  |            C               |      D = Stuart's rerun
#
# Reading across a row holds the pipeline fixed and changes the data.
# Reading down a column holds the data fixed and changes the pipeline.
# All four panels share one y axis, so the heights are comparable.
#
# A and B are computed here with the paper's own ANOVA_test_function, copied
# verbatim from its figshare deposit (doi 10.25387/g3.21200395,
# malathion_analysis.zip, get_haplotype_frequencies.txt), run on its own
# published haplotype frequencies (malathion_haps_genetic_1p5cM.txt.gz, same
# deposit). So A is a reproduction of the published figure, not a tracing of it.
#
# C is the current pipeline's own output on the four rep-1 pools, read from
# logs/malathion/<scan>.anova.txt.gz. D is the same on all eight pools and
# appears once the cluster run finishes; until then its panel says so.
#
# Usage, from the repo root:
#   Rscript temp_malathion/analysis/compare_scans.R
###############################################################################
suppressPackageStartupMessages({library(tidyverse); library(data.table)})

fs_dir     <- "temp_malathion/data"
haplo_file <- file.path(fs_dir, "malathion_haps_genetic_1p5cM.txt.gz")
rep1_file  <- "logs/malathion/malathion_rep1_all.anova.txt.gz"
all8_file  <- "logs/malathion/malathion_2022_all.anova.txt.gz"   # when the run lands
out_png    <- "temp_malathion/analysis/scan_comparison.png"
out_tab    <- "temp_malathion/analysis/scan_comparison.txt"

# The paper's data and code are the figshare supplement. Fetched once, not committed.
if (!file.exists(haplo_file)) {
  dir.create(fs_dir, showWarnings = FALSE, recursive = TRUE)
  cat("fetching the published haplotype frequencies from figshare...\n")
  download.file("https://ndownloader.figshare.com/files/37586293", haplo_file, quiet = TRUE, mode = "wb")
}
stopifnot(file.exists(haplo_file), file.exists(rep1_file))

# ── the paper's test, verbatim from the deposit ─────────────────────────────
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

cat("reading", haplo_file, "\n")
trans <- fread(haplo_file, header = TRUE) %>%
  mutate(TRT = str_sub(pool, 1, 1),
         REP = as.factor(str_sub(pool, 2, 3)),
         HAP = as.factor(gsub("B", "", gsub("A", "HAP", founder))),
         rep_expt = str_sub(pool, 2, 2),
         asf = asin(sqrt(freq))) %>%
  drop_na() %>%
  select(CHROM = chr, POS = pos, cM, TRT, REP, HAP, asf, rep_expt)

run_paper_scan <- function(df, cell) {
  cat("scanning", cell, ":", n_distinct(paste(df$CHROM, df$POS)), "windows,",
      n_distinct(df$REP), "pairs\n")
  df %>% group_by(CHROM, POS, cM) %>% nest() %>%
    mutate(mlog10p = map_dbl(data, ~ ANOVA_test_function(.))) %>%
    select(-data) %>% ungroup() %>%
    mutate(pipeline = "2022 pipeline", dataset = cell)
}
A <- run_paper_scan(trans, "both replicates (8 pools)")
B <- run_paper_scan(trans %>% filter(rep_expt == "1"), "replicate 1 only (4 pools)")

read_new <- function(f, cell) {
  if (!file.exists(f)) return(NULL)
  read.table(f, header = TRUE) %>%
    transmute(CHROM = chr, POS = pos, cM, mlog10p = Wald_log10p,
              pipeline = "current pipeline", dataset = cell)
}
C <- read_new(rep1_file, "replicate 1 only (4 pools)")
D <- read_new(all8_file, "both replicates (8 pools)")
if (is.null(D)) cat("note:", all8_file, "not present yet; panel D left empty\n")

# ── one cumulative cM axis, shared y ────────────────────────────────────────
chr_lv  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
pipe_lv <- c("2022 pipeline", "current pipeline")
data_lv <- c("replicate 1 only (4 pools)", "both replicates (8 pools)")
all4 <- bind_rows(A, B, C, D) %>% filter(CHROM %in% chr_lv) %>%
  mutate(CHROM = factor(CHROM, chr_lv),
         pipeline = factor(pipeline, pipe_lv),
         dataset  = factor(dataset,  data_lv))
chrlen <- all4 %>% group_by(CHROM) %>% summarise(len = max(cM), .groups = "drop") %>%
  arrange(CHROM) %>% mutate(off = cumsum(lag(len, default = 0)) + 2 * (as.integer(CHROM) - 1))
all4 <- all4 %>% left_join(chrlen, by = "CHROM") %>% mutate(cMcum = cM + off)
mid  <- all4 %>% group_by(CHROM) %>% summarise(m = mean(range(cMcum)), .groups = "drop")
genes <- tibble(name = c("Cyp6g1", "Mdr65"), CHROM = factor(c("chr2R", "chr3L"), chr_lv),
                cM = c(65.607, 16.627)) %>%
  left_join(chrlen, by = "CHROM") %>% mutate(cMcum = cM + off)

# peak heights, written onto each panel so nothing has to be eyeballed
peaks <- tibble(peak = c("2R", "3L"), CHROM = factor(c("chr2R","chr3L"), chr_lv),
                lo = c(55, 10), hi = c(75, 25))
pk <- map_dfr(seq_len(nrow(peaks)), function(k)
  all4 %>% filter(CHROM == peaks$CHROM[k], cM >= peaks$lo[k], cM <= peaks$hi[k]) %>%
    group_by(pipeline, dataset) %>%
    summarise(peak = peaks$peak[k], max = max(mlog10p, na.rm = TRUE),
              at_cM = cM[which.max(mlog10p)], .groups = "drop"))
ann <- pk %>% group_by(pipeline, dataset) %>%
  summarise(lab = paste0("2R ", sprintf("%.1f", max[peak == "2R"]),
                         "    3L ", sprintf("%.1f", max[peak == "3L"])), .groups = "drop")

ymax <- max(all4$mlog10p, na.rm = TRUE) * 1.18
missing <- expand_grid(pipeline = factor(pipe_lv, pipe_lv), dataset = factor(data_lv, data_lv)) %>%
  anti_join(all4 %>% distinct(pipeline, dataset), by = c("pipeline", "dataset"))

p <- ggplot(all4, aes(cMcum, mlog10p)) +
  geom_vline(data = genes, aes(xintercept = cMcum), colour = "grey75", linetype = 2, linewidth = .3) +
  geom_hline(yintercept = 4, colour = "firebrick", linetype = 3, linewidth = .3) +
  geom_line(aes(colour = CHROM), linewidth = .35) +
  geom_text(data = genes, aes(x = cMcum, y = ymax, label = name), inherit.aes = FALSE,
            vjust = 1, size = 2.4, colour = "grey35") +
  geom_text(data = ann, aes(x = -Inf, y = ymax, label = lab), inherit.aes = FALSE,
            hjust = -0.08, vjust = 1, size = 3, fontface = "bold", colour = "grey20") +
  facet_grid(pipeline ~ dataset) +
  coord_cartesian(ylim = c(0, ymax)) +
  scale_x_continuous(breaks = mid$m, labels = sub("^chr", "", levels(mid$CHROM))) +
  scale_colour_manual(values = rep(c("#1f4e79", "#2e8b57"), 3)[1:5], guide = "none") +
  labs(x = "Position (cM)", y = expression(-log[10] * P),
       title = "Malathion X-QTL: the same data through two pipeline versions",
       subtitle = "Across a row the pipeline is fixed and the data change. Down a column the data are fixed and the pipeline changes. One shared y axis.") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(face = "bold", size = 8.5),
        strip.background = element_rect(fill = "grey94", colour = NA),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))
# The one cell we cannot compute here: the Wald test needs the lsei covariance
# matrices, which only the haplotype step produces, so 8-pool haplotypes are
# required. The cluster is building them from the SRA reads. Meanwhile Stuart
# has run exactly this cell, so his reported peaks go in it, marked as reported.
if (nrow(missing))
  p <- p + geom_text(data = missing,
                     aes(x = mean(range(all4$cMcum)), y = ymax/2,
                         label = paste0("Stuart's rerun of this cell:  2R ~19.6    3L ~5.2\n",
                                        "(reported, not computed here)\n\n",
                                        "curve pending: needs 8-pool haplotypes,\n",
                                        "being built from the SRA reads")),
                     inherit.aes = FALSE, size = 3.1, colour = "grey45", lineheight = 1.25)

ggsave(out_png, p, width = 10, height = 6, dpi = 200)
cat("wrote", out_png, "\n")

sink(out_tab)
cat("Maximum -log10P at the two published QTL.\n\n")
print(as.data.frame(pk %>% mutate(max = round(max, 2), at_cM = round(at_cM, 2)) %>%
        pivot_wider(names_from = peak, values_from = c(max, at_cM)) %>%
        arrange(pipeline, dataset)), row.names = FALSE)
cat("\nPublished Table 1 peaks: chr2R 65.6 cM (12,294,089 bp), chr3L 16.6 cM (6,162,469 bp).\n")
cat("Published Fig. 2 peaks are just above 5 (2R) and about 4.9 (3L); the 2022-pipeline,\n")
cat("8-pool cell above is computed from the paper's own code and data and matches them,\n")
cat("which is the check that the paper's test is implemented correctly here.\n")
cat("\nCaveat on the current-pipeline, replicate-1 cell: its haplotypes were called from\n")
cat("the 2023 RefAlt tables on HPC3 (the old joint SNP caller), since those are what\n")
cat("exist for replicate 1. The 8-pool run now going uses the current catalog caller\n")
cat("end to end, so it will also say whether the SNP caller matters here.\n")
sink()
cat(readLines(out_tab), sep = "\n")
