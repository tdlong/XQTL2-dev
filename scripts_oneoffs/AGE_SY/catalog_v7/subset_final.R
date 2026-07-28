# subset_final.R — the final callset is a SUBSET of the loose v6 counts, not a recount.
# Take the already-merged loose RefAlt (process/AGE_SY_v6/Calls/RefAlt.<chr>.txt, from
# merge_dedup, multiallelic already dropped) and keep only positions passing the chosen
# filters (maxaf 3% + snpgap 20), read from snp_pass.tsv.gz. Writes v7 RefAlt.<chr>.txt.
#   Rscript subset_final.R <v6 Calls dir> <snp_pass.tsv.gz> <out Calls dir> <snpgap>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
v6 <- a[1]; passf <- a[2]; out <- a[3]; gap <- as.numeric(a[4])
dir.create(out, showWarnings = FALSE, recursive = TRUE)
chrs <- c("chrX","chr2L","chr2R","chr3L","chr3R")
pass <- fread(cmd = paste("zcat", shQuote(passf)))          # CHROM POS dist_indel pass03 pass05 pass07

for (chr in chrs) {
  f <- file.path(v6, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(f)) { cat("skip", chr, "\n"); next }
  ra   <- fread(cmd = paste("tr -s '\t ' ' ' <", shQuote(f)), header = TRUE, sep = " ")
  keep <- pass[CHROM == chr & pass03 == 1 & dist_indel >= gap, POS]
  final <- ra[POS %in% keep]
  fwrite(final, file.path(out, paste0("RefAlt.", chr, ".txt")), sep = "\t")
  cat(sprintf("%s: %d SNPs (of %d loose)\n", chr, nrow(final), nrow(ra)))
}
