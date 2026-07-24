# query_delta_gt2.R — QUERY the already-built delta_gt2 tables (the >2 tail).
# No re-merge: just reads process/AGE_SY_v4/Calls/delta_gt2.<chr>.tsv.gz (columns
# chr,pos,sample,total(=N_v3),delta,ref_v3,alt_v3,ref_v4,alt_v4) and answers:
#   1. how much of the >2 tail is high-coverage (max depth >=1000 -> the -d1000 cap)
#   2. within the NON-cap tail (max depth <1000): delta size + which samples
#   3. does the count disagreement also move the FREQUENCY (the thing that matters)
#
# Usage: Rscript query_delta_gt2.R <v4dir(=.../Calls)>
suppressPackageStartupMessages(library(data.table))
v4dir <- commandArgs(trailingOnly = TRUE)[1]
D <- rbindlist(lapply(Sys.glob(file.path(v4dir, "delta_gt2.*.tsv.gz")), fread))
D[, N3 := ref_v3 + alt_v3][, N4 := ref_v4 + alt_v4][, maxN := pmax(N3, N4)]
D[, dfreq := abs(ref_v3 / N3 - ref_v4 / N4)]
n <- nrow(D)
cat("delta>2 rows:", n, "\n")

cat("\n=== split by coverage (is the tail the -d 1000 cap?) ===\n")
cat(sprintf("  max depth >= 1000 (cap bites): %12d  %6.2f%%\n", D[maxN>=1000,.N], 100*D[maxN>=1000,.N]/n))
cat(sprintf("  max depth <  1000 (no cap)   : %12d  %6.2f%%\n", D[maxN<1000,.N],  100*D[maxN<1000,.N]/n))

lo <- D[maxN < 1000]
cat(sprintf("\n=== the NON-cap tail (%d rows) — delta size ===\n", nrow(lo)))
lo[, .N, by = .(delta_bin = cut(delta, c(2,5,10,20,Inf), labels=c("3-5","6-10","11-20",">20")))][order(delta_bin)][, cat(sprintf("  %-6s %12d  %6.2f%%\n", delta_bin, N, 100*N/nrow(lo))), by=.(delta_bin,N)]

cat("\n=== NON-cap tail: top 15 samples by count ===\n")
print(lo[, .N, by = sample][order(-N)][1:15])

cat("\n=== does the count disagreement also move FREQUENCY? (all >2 rows) ===\n")
D[, .N, by = .(dfreq_bin = cut(dfreq, c(-1,0.01,0.02,0.05,0.10,1), labels=c("<0.01","0.01-0.02","0.02-0.05","0.05-0.10",">0.10")))][order(dfreq_bin)][, cat(sprintf("  dfreq %-10s %12d  %6.2f%%\n", dfreq_bin, N, 100*N/n)), by=.(dfreq_bin,N)]
cat("\n(dfreq = |freq_v3 - freq_v4|; big count delta but tiny dfreq = harmless, ratio preserved)\n")
