# agreement_delta.R — per-(SNP,sample) count agreement, v3 vs v4 (same BAM).
#
#   total = REF_v3 + ALT_v3
#   delta = |REF_v3 - REF_v4| + |ALT_v3 - ALT_v4|   (reads that moved; 0 = exact)
#
# delta 0/1/2 are the expected nudge from the callers' different filters
# (v3: BAQ-on, -q 0, default -Q, -d 1000; v4: BAQ-off, -q 20, -Q 20, --max-depth
# 2000) — we only TALLY those. delta > 2 is the tail worth inspecting — we KEEP it.
#
# OUTPUTS
#  (1) log: proportion of SNP*samples at delta 0, 1, 2, >2 (overall + per chr),
#      and a summary of the >2 tail.
#  (2) <v4dir>/delta_gt2.<chr>.tsv.gz: chr,pos,sample,total,delta,ref_v3,alt_v3,
#      ref_v4,alt_v4 for every pair with delta > 2 (the reusable table to query).
#
# Usage: Rscript agreement_delta.R <v3dir> <v4dir(=.../Calls)>
suppressPackageStartupMessages(library(data.table))
args  <- commandArgs(trailingOnly = TRUE)
v3dir <- args[1]; v4dir <- args[2]
chrs  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
readRA <- function(f) fread(cmd = paste("tr -s '\\t ' ' ' <", shQuote(f)), header = TRUE, sep = " ")

c0 <- 0; c1 <- 0; c2 <- 0; cg <- 0        # counts of delta 0 / 1 / 2 / >2
gt <- integer(4); names(gt) <- c("3-5","6-10","11-20",">20")  # within the >2 tail
top <- data.table(); perchr <- list()

for (chr in chrs) {
  fa <- file.path(v3dir, paste0("RefAlt.", chr, ".txt"))
  fb <- file.path(v4dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(fa) || !file.exists(fb)) { cat("skip", chr, "\n"); next }
  m <- merge(readRA(fa), readRA(fb), by = c("CHROM", "POS"), suffixes = c(".v3", ".v4"))
  samps <- sub("\\.v4$", "", sub("^REF_", "", grep("^REF_.*\\.v4$", names(m), value = TRUE)))

  keep <- vector("list", length(samps)); k0<-0; k1<-0; k2<-0; kg<-0
  for (i in seq_along(samps)) { s <- samps[i]
    r3 <- m[[paste0("REF_",s,".v3")]]; a3 <- m[[paste0("ALT_",s,".v3")]]
    r4 <- m[[paste0("REF_",s,".v4")]]; a4 <- m[[paste0("ALT_",s,".v4")]]
    if (is.null(r3) || is.null(r4)) next
    d <- abs(r3 - r4) + abs(a3 - a4); N <- r3 + a3
    k0 <- k0 + sum(d==0); k1 <- k1 + sum(d==1); k2 <- k2 + sum(d==2)
    g <- d > 2; kg <- kg + sum(g)
    if (any(g)) keep[[i]] <- data.table(chr=chr, pos=m$POS[g], sample=s, total=N[g], delta=d[g],
        ref_v3=r3[g], alt_v3=a3[g], ref_v4=r4[g], alt_v4=a4[g])
  }
  c0<-c0+k0; c1<-c1+k1; c2<-c2+k2; cg<-cg+kg
  perchr[[chr]] <- c(pairs=k0+k1+k2+kg, d0=k0, d1=k1, d2=k2, dgt2=kg)
  G <- rbindlist(keep)
  if (nrow(G)) {
    gt["3-5"]<-gt["3-5"]+sum(G$delta<=5); gt["6-10"]<-gt["6-10"]+sum(G$delta>=6&G$delta<=10)
    gt["11-20"]<-gt["11-20"]+sum(G$delta>=11&G$delta<=20); gt[">20"]<-gt[">20"]+sum(G$delta>20)
    fwrite(G, file.path(v4dir, paste0("delta_gt2.", chr, ".tsv.gz")), sep = "\t")
    top <- rbindlist(list(top, G[order(-delta)][1:min(200,.N)]))[order(-delta)][1:min(200,.N)]
  }
  rm(m, G, keep); gc()
}

tot <- c0 + c1 + c2 + cg
p <- function(x) sprintf("%14d  %6.2f%%", x, 100*x/tot)
cat("=== per-(SNP,sample) count agreement, v3 vs v4 (shared SNPs & samples) ===\n")
cat("  delta == 0 (exact) :", p(c0), "\n")
cat("  delta == 1         :", p(c1), "\n")
cat("  delta == 2         :", p(c2), "\n")
cat("  delta  > 2  (KEPT) :", p(cg), "\n")
cat("  total pairs        :", sprintf("%14d", tot), "\n")

cat("\n=== per chromosome (pairs / d=0 / d=1 / d=2 / d>2) ===\n")
cat(sprintf("%-8s %12s %12s %12s %12s %12s\n","chr","pairs","d0","d1","d2","d>2"))
for (chr in names(perchr)) { v<-perchr[[chr]]
  cat(sprintf("%-8s %12d %12d %12d %12d %12d\n", chr, v["pairs"],v["d0"],v["d1"],v["d2"],v["dgt2"])) }

cat("\n=== the delta>2 tail, by size ===\n")
for (nm in names(gt)) cat(sprintf("  delta %-6s %12d  %6.2f%% of the >2 tail\n", nm, gt[nm], if(cg>0) 100*gt[nm]/cg else 0))
cat("\n=== 20 largest deltas (what big discrepancies look like) ===\n")
if (nrow(top)) print(head(top[order(-delta)], 20))
cat(sprintf("\nkept table: %s/delta_gt2.<chr>.tsv.gz  (query this, no rebuild)\n", v4dir))
