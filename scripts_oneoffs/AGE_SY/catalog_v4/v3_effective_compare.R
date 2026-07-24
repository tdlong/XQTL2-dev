# v3_effective_compare.R — three-way SNP-set comparison:
#   v3_raw       every SNP in v3 RefAlt (QUAL>59)
#   v3_effective v3 after the good_SNPs founder filter (what the haplotyper sees)
#   catalog      the v4 founder-catalog SNP set
#
# good_SNPs is replicated FAITHFULLY from REFALT2haps.code.R:169-176 (founders
# only, per SNP): zeros==0 (every founder covered) & notfixed==0 (every covered
# founder near-fixed at 0.03/0.97). The 'informative' term there is an OR that is
# true for every kept SNP, so it is a no-op and omitted.
#
# Usage:
#   Rscript v3_effective_compare.R <v3dir> <catalog.tsv.gz> <parfile.R>
suppressPackageStartupMessages(library(data.table))
args    <- commandArgs(trailingOnly = TRUE)
v3dir   <- args[1]
catgz   <- args[2]
parfile <- args[3]
chrs    <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

source(parfile)                      # defines `founders`
cat("founders:", paste(founders, collapse = " "), "\n\n")

catdt <- fread(cmd = paste("zcat", shQuote(catgz)), header = FALSE)  # CHROM POS REF,ALT
setnames(catdt, c("CHROM", "POS", "REFALT"))

cat(sprintf("%-7s %10s %10s %10s | %9s %9s %9s\n",
            "chr", "v3_raw", "v3_eff", "catalog", "cat&eff", "cat_only", "eff_only"))
tot <- c(raw = 0, eff = 0, cat = 0, inter = 0, conly = 0, eonly = 0)
for (chr in chrs) {
  f <- file.path(v3dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(f)) { cat(sprintf("%-7s (missing %s)\n", chr, f)); next }
  dt <- fread(f, header = TRUE)      # CHROM POS REF_<name> ALT_<name> ...

  zeros <- integer(nrow(dt)); notfixed <- integer(nrow(dt))
  for (fd in founders) {
    rc <- dt[[paste0("REF_", fd)]]; ac <- dt[[paste0("ALT_", fd)]]
    if (is.null(rc) || is.null(ac)) stop(paste("missing founder column for", fd, "in", f))
    N <- rc + ac; freq <- rc / N
    zeros    <- zeros    + as.integer(N == 0)
    nf <- as.integer(N != 0 & freq > 0.03 & freq < 0.97); nf[is.na(nf)] <- 0L
    notfixed <- notfixed + nf
  }
  good <- zeros == 0 & notfixed == 0

  effset <- dt[good, POS]
  catset <- catdt[CHROM == chr, POS]
  inter  <- length(intersect(catset, effset))
  conly  <- length(setdiff(catset, effset))
  eonly  <- length(setdiff(effset, catset))
  cat(sprintf("%-7s %10d %10d %10d | %9d %9d %9d\n",
              chr, nrow(dt), length(effset), length(catset), inter, conly, eonly))
  tot <- tot + c(nrow(dt), length(effset), length(catset), inter, conly, eonly)
}
cat(sprintf("%-7s %10d %10d %10d | %9d %9d %9d\n",
            "TOTAL", tot["raw"], tot["eff"], tot["cat"], tot["inter"], tot["conly"], tot["eonly"]))
cat("\nread: cat_only = catalog SNPs NOT in v3-effective; eff_only = good_SNPs the catalog dropped.\n")
cat("prediction: catalog ~subset of v3_eff (cat_only small); eff_only = founder-monomorphic SNPs good_SNPs keeps but catalog drops as non-segregating.\n")
