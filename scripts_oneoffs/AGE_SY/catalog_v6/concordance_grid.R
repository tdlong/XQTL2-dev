# concordance_grid.R — STEP 4. For each of the 12 filter combos (indel 0/10/20/50 x
# maxaf 3/5/7%), what fraction of (SNP,sample) pairs differ from v3 by >2% and >5% in
# ALT frequency. Reuses the RefAlt merge from agreement_delta.R.
#
#   Rscript concordance_grid.R <v6 Calls dir> <v3 dir> <snp_pass.tsv.gz>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); v6dir <- a[1]; v3dir <- a[2]; passf <- a[3]
chrs <- c("chrX","chr2L","chr2R","chr3L","chr3R")
readRA <- function(f) fread(cmd = paste("tr -s '\\t ' ' ' <", shQuote(f)), header = TRUE, sep = " ")

pass <- fread(cmd = paste("zcat", shQuote(passf)))            # CHROM POS dist_indel pass03 pass05 pass07
G <- c(0,10,20,50); M <- c("pass03","pass05","pass07"); Mlab <- c("3%","5%","7%")

# accumulators: for each (indel g, maxaf m): n pairs, n>2%, n>5%
n  <- array(0, c(length(G), length(M)))
g2 <- array(0, c(length(G), length(M)))
g5 <- array(0, c(length(G), length(M)))

for (chr in chrs) {
  f6 <- file.path(v6dir, paste0("RefAlt.", chr, ".txt"))
  f3 <- file.path(v3dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(f6) || !file.exists(f3)) { cat("skip", chr, "\n"); next }
  m <- merge(readRA(f3), readRA(f6), by = c("CHROM","POS"), suffixes = c(".v3",".v6"))
  m <- merge(m, pass[CHROM == chr], by = c("CHROM","POS"))     # attach dist_indel + pass flags
  samps <- sub("\\.v6$","", sub("^REF_","", grep("^REF_.*\\.v6$", names(m), value = TRUE)))
  for (s in samps) {
    r3 <- m[[paste0("REF_",s,".v3")]]; a3 <- m[[paste0("ALT_",s,".v3")]]
    r6 <- m[[paste0("REF_",s,".v6")]]; a6 <- m[[paste0("ALT_",s,".v6")]]
    if (is.null(r3) || is.null(r6)) next
    N3 <- r3 + a3; N6 <- r6 + a6
    ok <- N3 > 0 & N6 > 0                                      # both callers have reads
    df <- abs(a6/N6 - a3/N3)
    for (im in seq_along(M)) {
      pm <- m[[M[im]]] == 1
      for (ig in seq_along(G)) {
        sel <- ok & pm & (m$dist_indel >= G[ig])
        ns <- sum(sel)
        n[ig,im]  <- n[ig,im]  + ns
        g2[ig,im] <- g2[ig,im] + sum(df[sel] > 0.02)
        g5[ig,im] <- g5[ig,im] + sum(df[sel] > 0.05)
      }
    }
  }
  rm(m); gc()
}

cat("=== CONCORDANCE v6(BAQ-on catalog) vs v3, per filter combo ===\n")
cat("fraction of (SNP,sample) pairs whose ALT freq differs from v3 by >2%  /  >5%\n\n")
for (im in seq_along(M)) {
  cat(sprintf("--- maxaf %s ---\n", Mlab[im]))
  cat(sprintf("  %-8s %12s %12s %12s\n", "indel>=", "pairs", ">2%", ">5%"))
  for (ig in seq_along(G)) {
    cat(sprintf("  %-8d %12d %11.2f%% %11.2f%%\n",
        G[ig], n[ig,im], 100*g2[ig,im]/max(1,n[ig,im]), 100*g5[ig,im]/max(1,n[ig,im])))
  }
  cat("\n")
}
