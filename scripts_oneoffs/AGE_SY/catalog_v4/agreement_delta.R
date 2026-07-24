# agreement_delta.R — per-(SNP,sample) count agreement, v3 vs v4 (same BAM).
#   N3 = REF_v3+ALT_v3   N4 = REF_v4+ALT_v4
#   delta = |REF_v3-REF_v4| + |ALT_v3-ALT_v4|   (0 = exact)
# Tallies delta 0/1/2/>2, STRATIFIED by max(N3,N4): <1000 (neither caller hit its
# depth cap; isolates -q/-Q/BAQ) vs >=1000 (v3's -d 1000 vs v4's 2000 can bite).
# Keeps the >2 rows to delta_gt2.<chr>.tsv.gz (chr,pos,sample,N3,N4,delta,ref/alt).
#
# Usage: Rscript agreement_delta.R <v3dir> <v4dir(=.../Calls)>
suppressPackageStartupMessages(library(data.table))
args  <- commandArgs(trailingOnly = TRUE); v3dir <- args[1]; v4dir <- args[2]
chrs  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
readRA <- function(f) fread(cmd = paste("tr -s '\\t ' ' ' <", shQuote(f)), header = TRUE, sep = " ")

lo <- integer(4); hi <- integer(4)          # counts of delta 0,1,2,>2 per depth stratum
top <- data.table()
for (chr in chrs) {
  fa <- file.path(v3dir, paste0("RefAlt.", chr, ".txt")); fb <- file.path(v4dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(fa) || !file.exists(fb)) { cat("skip", chr, "\n"); next }
  m <- merge(readRA(fa), readRA(fb), by = c("CHROM", "POS"), suffixes = c(".v3", ".v4"))
  samps <- sub("\\.v4$", "", sub("^REF_", "", grep("^REF_.*\\.v4$", names(m), value = TRUE)))
  keep <- vector("list", length(samps))
  for (i in seq_along(samps)) { s <- samps[i]
    r3 <- m[[paste0("REF_",s,".v3")]]; a3 <- m[[paste0("ALT_",s,".v3")]]
    r4 <- m[[paste0("REF_",s,".v4")]]; a4 <- m[[paste0("ALT_",s,".v4")]]
    if (is.null(r3) || is.null(r4)) next
    d <- abs(r3-r4) + abs(a3-a4); N3 <- r3+a3; N4 <- r4+a4
    b <- pmin(d, 3L) + 1L                   # 1..4 = delta 0,1,2,>2
    ishi <- pmax(N3, N4) >= 1000
    lo <- lo + tabulate(b[!ishi], 4L); hi <- hi + tabulate(b[ishi], 4L)
    g <- d > 2
    if (any(g)) keep[[i]] <- data.table(chr=chr, pos=m$POS[g], sample=s, N3=N3[g], N4=N4[g],
        delta=d[g], ref_v3=r3[g], alt_v3=a3[g], ref_v4=r4[g], alt_v4=a4[g])
  }
  G <- rbindlist(keep)
  if (nrow(G)) { fwrite(G, file.path(v4dir, paste0("delta_gt2.", chr, ".tsv.gz")), sep = "\t")
    top <- rbindlist(list(top, G[order(-delta)][1:min(50,.N)]))[order(-delta)][1:min(50,.N)] }
  rm(m, G, keep); gc()
}

show <- function(v, lab) { t <- sum(v)
  cat(sprintf("--- %s : %d pairs ---\n", lab, t))
  for (i in 1:4) cat(sprintf("  delta %-3s %14d  %6.2f%%\n", c("0","1","2",">2")[i], v[i], 100*v[i]/t)) }
cat("=== per-(SNP,sample) count agreement, v3 vs v4 ===\n")
show(lo, "max(depth) < 1000  (no -d cap effect; pure -q/-Q/BAQ)")
show(hi, "max(depth) >= 1000 (v3 -d1000 vs v4 2000 can bite)")
show(lo + hi, "ALL")
cat("\n=== 20 largest deltas ===\n"); if (nrow(top)) print(head(top[order(-delta)], 20))
cat(sprintf("\nkept: %s/delta_gt2.<chr>.tsv.gz\n", v4dir))
