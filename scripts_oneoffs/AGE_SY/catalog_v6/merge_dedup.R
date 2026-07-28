# merge_dedup.R — workaround for XQTL2 #24: build RefAlt.<chr>.txt from the per-sample
# catalog counts, DROPPING multiallelic positions (same POS >1 row) that make
# catalog_merge's (CHROM,POS) join explode. The counting already ran; this only merges.
#   Rscript merge_dedup.R <calls dir>   (reads <calls dir>/counts/*.tsv.gz)
suppressPackageStartupMessages(library(data.table))
callsdir <- commandArgs(trailingOnly = TRUE)[1]
chrs  <- c("chrX","chr2L","chr2R","chr3L","chr3R")
files <- list.files(file.path(callsdir, "counts"), pattern = "\\.tsv\\.gz$", full.names = TRUE)
stopifnot(length(files) > 0)

for (chr in chrs) {
  merged <- NULL; ndrop <- 0L
  for (f in files) {
    dt <- fread(cmd = paste("gzip -dc", shQuote(f)), sep = "\t", header = TRUE)
    dt <- dt[CHROM == chr]
    dup <- dt[, .N, by = POS][N > 1, POS]          # multiallelic positions in this catalog
    if (length(dup)) { dt <- dt[!POS %in% dup]; ndrop <- length(dup) }
    merged <- if (is.null(merged)) dt else merge(merged, dt, by = c("CHROM","POS"), all = TRUE)
  }
  if (is.null(merged) || !nrow(merged)) { cat("skip", chr, "\n"); next }
  cc <- setdiff(names(merged), c("CHROM","POS"))
  merged[, (cc) := lapply(.SD, function(x) fifelse(is.na(x), 0L, as.integer(x))), .SDcols = cc]
  setorder(merged, POS)
  fwrite(merged, file.path(callsdir, paste0("RefAlt.", chr, ".txt")), sep = "\t")
  cat(sprintf("%s: %d SNPs written, %d multiallelic positions dropped\n", chr, nrow(merged), ndrop))
}
