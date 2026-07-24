# v3_effective_compare.R — three-way SNP-set comparison + breakdown of both tails.
#
#   v3_raw       every SNP in v3 RefAlt (QUAL>59)
#   v3_effective v3 after the good_SNPs founder filter (what the haplotyper sees)
#   catalog      the v4 founder-catalog SNP set
#
# good_SNPs replicated from REFALT2haps.code.R:169-176 (founders only, per SNP):
# zeros==0 (every founder covered) & notfixed==0 (every covered founder near-fixed
# at 0.03/0.97). 'informative' there is an always-true OR -> no-op, omitted.
#
# Then breaks down the two disagreement tails:
#   cat_only  = catalog SNPs NOT in v3-effective, split into
#               QUAL_missed  (not in v3_raw at all -> the catalog's real gain)
#               foundr_disag (in v3_raw but good_SNPs dropped -> pooled vs direct founder call)
#   eff_only  = good_SNPs SNPs the catalog dropped, split into
#               monomorphic  (all founders fixed for the SAME allele -> catalog's segregating rule; uninformative anyway)
#               segregating  (>=1 REF-fixed & >=1 ALT-fixed founder, yet dropped -> near-indel snpgap / measurement; the real question)
#
# Usage:
#   Rscript v3_effective_compare.R <v3dir> <catalog.tsv.gz> <parfile.R>
suppressPackageStartupMessages(library(data.table))
args    <- commandArgs(trailingOnly = TRUE)
v3dir   <- args[1]; catgz <- args[2]; parfile <- args[3]
chrs    <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
MAXAF   <- 0.03

source(parfile)                      # defines `founders`
cat("founders:", paste(founders, collapse = " "), "\n\n")

catdt <- fread(cmd = paste("zcat", shQuote(catgz)), header = FALSE)  # CHROM POS REF,ALT
setnames(catdt, c("CHROM", "POS", "REFALT"))

A <- list(); B <- list(); C <- list()   # per-chr rows for the three tables
for (chr in chrs) {
  f <- file.path(v3dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(f)) { cat(sprintf("skip %s (missing %s)\n", chr, f)); next }
  # v3 RefAlt: TAB header, SPACE data (bam2bcf2REFALT quirk) -> normalize to spaces.
  dt <- fread(cmd = paste("tr '\\t' ' ' <", shQuote(f)), header = TRUE, sep = " ")

  zeros <- integer(nrow(dt)); notfixed <- integer(nrow(dt))
  nreff <- integer(nrow(dt)); naltf <- integer(nrow(dt))   # founders fixed ref / alt
  for (fd in founders) {
    rc <- dt[[paste0("REF_", fd)]]; ac <- dt[[paste0("ALT_", fd)]]
    if (is.null(rc) || is.null(ac)) stop(paste("missing founder column for", fd, "in", f))
    N <- rc + ac; freq <- rc / N                           # freq = REF fraction
    zeros    <- zeros    + as.integer(N == 0)
    nf <- as.integer(N != 0 & freq > MAXAF & freq < 1 - MAXAF); nf[is.na(nf)] <- 0L
    notfixed <- notfixed + nf
    rf <- as.integer(N != 0 & freq >= 1 - MAXAF); rf[is.na(rf)] <- 0L
    af <- as.integer(N != 0 & freq <= MAXAF);     af[is.na(af)] <- 0L
    nreff <- nreff + rf; naltf <- naltf + af
  }
  good        <- zeros == 0 & notfixed == 0                 # v3-effective mask
  segregating <- nreff >= 1 & naltf >= 1                    # >=1 ref-fixed AND >=1 alt-fixed founder

  rawPOS <- dt$POS
  effPOS <- dt$POS[good]
  catPOS <- catdt[CHROM == chr, POS]
  inCat  <- effPOS %in% catPOS                              # for eff_only we need the NOT-in-catalog ones
  effOnlyPOS <- effPOS[!inCat]

  # A. three-way
  inter <- length(intersect(catPOS, effPOS))
  catOnly <- setdiff(catPOS, effPOS)
  A[[chr]] <- c(raw = length(rawPOS), eff = length(effPOS), cat = length(catPOS),
                inter = inter, cat_only = length(catOnly), eff_only = length(effOnlyPOS))

  # B. cat_only split: not in v3_raw (QUAL missed) vs in raw but good_SNPs dropped
  qmiss <- length(setdiff(catOnly, rawPOS))
  B[[chr]] <- c(cat_only = length(catOnly), QUAL_missed = qmiss, foundr_disag = length(catOnly) - qmiss)

  # C. eff_only split: monomorphic (non-segregating) vs still-segregating
  segByPOS <- segregating[good][!inCat]                    # segregation status of the eff_only SNPs
  C[[chr]] <- c(eff_only = length(effOnlyPOS),
                monomorphic = sum(!segByPOS), segregating = sum(segByPOS))
}

prtab <- function(L, hdr, cols) {
  cat("\n=== ", hdr, " ===\n", sep = "")
  row <- function(lab, vals) cat(formatC(lab, width = 8, flag = "-"),
                                 paste(formatC(vals, width = 12), collapse = ""), "\n")
  row("chr", cols)
  tot <- NULL
  for (chr in names(L)) { v <- L[[chr]]; tot <- if (is.null(tot)) v else tot + v; row(chr, v) }
  row("TOTAL", tot)
}

prtab(A, "A. three-way SNP sets", c("v3_raw","v3_eff","catalog","cat&eff","cat_only","eff_only"))
prtab(B, "B. cat_only breakdown (catalog SNPs not in v3-effective)", c("cat_only","QUAL_missed","foundr_disag"))
prtab(C, "C. eff_only breakdown (good_SNPs the catalog dropped)",    c("eff_only","monomorphic","segregating"))
cat("\nread:\n")
cat("  B QUAL_missed = catalog SNPs QUAL>59 never called (the gain); foundr_disag = pooled(v3) vs direct(catalog) founder-call disagreement.\n")
cat("  C monomorphic = all founders fixed same allele (catalog drops as non-segregating; useless for haps anyway).\n")
cat("  C segregating = real founder-segregating SNPs the catalog still dropped -> near-indel (snpgap) or measurement; the ones to inspect.\n")
