# coverage.R — observed depth at catalog SNPs, over the 10 replicates kept.
#
#   Rscript temp_aging/coverage.R
#
# RUN ON HPC3, from /dfs7/adl/tdlong/fly_pool/XQTL2-dev -- process/ lives there.
#
#   Rscript pipeline/scripts/refalt_qc.R \
#       --dir     process/AGE_SY \
#       --parfile helpfiles/AGE_SY/AGE_SY_haplotype_parameters_size75k.R
#   Rscript temp_aging/coverage.R
#
# The first writes process/AGE_SY/Calls/refalt_qc.txt (one row per sample per
# chromosome, median REF+ALT over catalog sites); this reads it. Neither needs a
# SLURM job -- refalt_qc.R reads the RefAlt tables once, coverage.R reads one
# small table. To keep the output with the rest, run it through run_numbers.sh
# instead and scp back numbers/coverage.txt, which is a few hundred bytes:
#
#   bash temp_aging/run_numbers.sh          # on HPC3
#   scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/temp_aging/numbers/coverage.txt \
#       temp_aging/numbers/
#
# The catalog and the counts were built on all 12 replicates. The scans use 10 --
# 8 and 9 dropped as the May 2023 cage. This reports the 10 that are kept: 60
# libraries, 3 groups x 2 sexes x 10 cages, controls included, founders excluded.
#
# chrX is reported apart from the autosomes and not averaged into them. Males are
# hemizygous, so male chrX depth is about half male autosomal depth; a single
# genome-wide mean would hide that, and paragraph 3 turns on an X-vs-autosome
# contrast.

suppressMessages(library(tidyverse))

QC <- "process/AGE_SY/Calls/refalt_qc.txt"
if (!file.exists(QC)) stop(
  "missing ", QC, "\nRun pipeline/scripts/refalt_qc.R on HPC3 and scp it back; ",
  "see the header of this script.", call. = FALSE)

DROP_REP <- c(8, 9)

d <- read.table(QC, header = TRUE, sep = "\t") %>% as_tibble() %>%
  # sample names are <group>_R<rep>_<sex>; founders (B1..B7, AB8) do not match
  filter(str_detect(sample, "^(AgeSY10|AgeSY20|Con)_R[0-9]+_[FM]$")) %>%
  mutate(group = str_extract(sample, "^[^_]+"),
         rep   = as.integer(str_extract(sample, "(?<=_R)[0-9]+")),
         sex   = str_extract(sample, "[FM]$"),
         arm   = if_else(chr == "chrX", "chrX", "autosome")) %>%
  filter(!rep %in% DROP_REP)

stopifnot(n_distinct(d$sample) == 60, n_distinct(d$rep) == 10)

# One depth per sample per arm: median over that sample's chromosomes, so the
# four autosomal arms are not counted as four samples.
per_sample <- d %>% group_by(sample, group, sex, arm) %>%
  summarise(depth = median(median_depth), .groups = "drop")

cat("libraries: ", n_distinct(d$sample), " (", n_distinct(d$rep),
    " replicates x 3 groups x 2 sexes)\n\n", sep = "")

cat("=== median depth at catalog SNPs, across the 60 libraries ===\n")
per_sample %>% group_by(arm, sex) %>%
  summarise(n = n(), mean = round(mean(depth), 1), median = round(median(depth), 1),
            min = round(min(depth), 1), max = round(max(depth), 1), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n=== by group ===\n")
per_sample %>% group_by(arm, sex, group) %>%
  summarise(n = n(), mean = round(mean(depth), 1),
            min = round(min(depth), 1), max = round(max(depth), 1), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n=== the number for the methods ===\n")
aut <- per_sample %>% filter(arm == "autosome")
cat(sprintf("autosomal median depth per library: mean %.0fx, range %.0f-%.0fx (n=%d)\n",
            mean(aut$depth), min(aut$depth), max(aut$depth), nrow(aut)))
xm <- per_sample %>% filter(arm == "chrX", sex == "M")
xf <- per_sample %>% filter(arm == "chrX", sex == "F")
cat(sprintf("chrX: females mean %.0fx, males mean %.0fx (%.2f of female)\n",
            mean(xf$depth), mean(xm$depth), mean(xm$depth) / mean(xf$depth)))

# How many libraries sit below a given depth, for the methods sentence. Reported
# for the autosomes; male chrX is listed separately because hemizygosity puts it
# near 72x by construction, so counting it against the same threshold would read
# as a library problem when it is the expected halving.
cat("\n=== libraries below depth thresholds ===\n")
for (t in c(75, 100)) {
  cat(sprintf("  autosomal median < %3dx : %2d of %d libraries\n",
              t, sum(aut$depth < t), nrow(aut)))
}
for (t in c(75, 100)) {
  cat(sprintf("  male chrX      < %3dx : %2d of %d male libraries\n",
              t, sum(xm$depth < t), nrow(xm)))
}

cat("\n=== thin libraries (worst 5 on autosomal depth) ===\n")
aut %>% arrange(depth) %>% head(5) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\n=== anything the pipeline flagged ===\n")
f <- d %>% filter(flag != "OK") %>% select(sample, chr, median_depth, pct_zero,
                                           pct_window_covered, flag)
if (nrow(f) == 0) cat("none of the 60 flagged on any chromosome.\n") else
  print(as.data.frame(f), row.names = FALSE)
