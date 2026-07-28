# compare_summary.R — a better v3-vs-v6 summary than a single mean:
#   (1) SNP-set overlap as plain counts (v3-only / shared / v6-only), per chr + total.
#   (2) the JOINT distribution of |ALT-freq difference| vs COVERAGE at shared (SNP,sample)
#       -- binned by depth, with the spread (median/p90/p99, fraction over thresholds) and
#       the binomial-expected |diff| per bin, so we can see if agreement tracks sampling.
#   Also writes a plot (|freq diff| percentiles vs coverage) to figures/.
#
#   Rscript compare_summary.R <v3 dir> <v6 Calls dir>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE); v3 <- a[1]; v6 <- a[2]
chrs <- c("chrX","chr2L","chr2R","chr3L","chr3R")
readRA <- function(f) fread(cmd = paste("tr -s '\t ' ' ' <", shQuote(f)), header = TRUE, sep = " ")

ov <- data.table(); D <- vector("list", length(chrs))
for (ci in seq_along(chrs)) { chr <- chrs[ci]
  fa <- file.path(v3, paste0("RefAlt.", chr, ".txt")); fb <- file.path(v6, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(fa) || !file.exists(fb)) { cat("skip", chr, "\n"); next }
  A <- readRA(fa); B <- readRA(fb)
  posA <- unique(A[, .(CHROM, POS)]); posB <- unique(B[, .(CHROM, POS)])
  sh <- fintersect(posA, posB)
  ov <- rbind(ov, data.table(chr = chr, v3_only = nrow(posA) - nrow(sh),
                             shared = nrow(sh), v6_only = nrow(posB) - nrow(sh)))
  m <- merge(A, B, by = c("CHROM","POS"), suffixes = c(".v3",".v6"))
  samps <- sub("\\.v6$","", sub("^REF_","", grep("^REF_.*\\.v6$", names(m), value = TRUE)))
  L <- vector("list", length(samps))
  for (i in seq_along(samps)) { s <- samps[i]
    r3<-m[[paste0("REF_",s,".v3")]]; a3<-m[[paste0("ALT_",s,".v3")]]
    r6<-m[[paste0("REF_",s,".v6")]]; a6<-m[[paste0("ALT_",s,".v6")]]
    if (is.null(r3)||is.null(r6)) next
    N3<-r3+a3; N6<-r6+a6; ok<-N3>0 & N6>0
    f3<-(a3/N3)[ok]; f6<-(a6/N6)[ok]; pm<-(f3+f6)/2
    L[[i]] <- data.table(cov = pmin(N3,N6)[ok], d = abs(f6-f3), sd = f6-f3,
                         p = f6, maf = pmin(pm, 1-pm))
  }
  D[[ci]] <- rbindlist(L); rm(m); gc()
}
DT <- rbindlist(D)

cat("=== SNP-set overlap (positions) ===\n")
ovt <- rbind(ov, data.table(chr="TOTAL", v3_only=sum(ov$v3_only), shared=sum(ov$shared), v6_only=sum(ov$v6_only)))
print(ovt)

cat("\n=== |ALT-freq diff| vs coverage at shared (SNP,sample) ===\n")
brk <- c(0,10,25,50,100,200,400,Inf)
DT[, cbin := cut(cov, brk, right = FALSE,
   labels=c("1-9","10-24","25-49","50-99","100-199","200-399","400+"))]
tab <- DT[, .(n=.N, median=round(median(d),4), p90=round(quantile(d,.9),4),
              p99=round(quantile(d,.99),4),
              pct_gt2=round(100*mean(d>0.02),2), pct_gt5=round(100*mean(d>0.05),2),
              binom_exp=round(mean(sqrt(2*p*(1-p)/cov)),4)), by = cbin][order(cbin)]
print(tab)
cat("\nbinom_exp = mean sqrt(2p(1-p)/cov) — the |diff| expected from sampling alone.\n")
cat("If median tracks binom_exp, the two callers agree within counting noise.\n")

# signed diff by minor-allele-freq bin (is any freq range biased? esp. rare alleles / #22)
DT[, mbin := cut(maf, c(0,.05,.1,.2,.35,.5), include.lowest=TRUE,
   labels=c("MAF 0-.05",".05-.1",".1-.2",".2-.35",".35-.5"))]
cat("\n=== signed diff (freq_v6 - freq_v3) by minor-allele-freq bin ===\n")
print(DT[, .(n=.N, mean_sd=round(mean(sd),4), median_sd=round(median(sd),4),
             p05=round(quantile(sd,.05),4), p95=round(quantile(sd,.95),4),
             pct_gt2=round(100*mean(d>.02),2)), by=mbin][order(mbin)])

# write to logs/figures/ — logs/ already syncs, so the plots come back with cluster_sync
FIG <- "logs/figures"; dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(FIG,"compare_v3_v6_freqdiff_vs_cov.png"), width=7, height=5, dpi=120,
  ggplot(DT[sample(.N, min(.N,3e5))], aes(cov, d)) + geom_hex(bins=60) +
    scale_x_log10() + scale_y_sqrt() + theme_minimal() +
    labs(x="coverage min(N3,N6)", y="|freq_v6 - freq_v3|", title="freq diff vs coverage"))

# THE density-by-frequency-bin plot: 5 normalized curves, common bin can't dominate
ggsave(file.path(FIG,"compare_v3_v6_freqdiff_density_by_maf.png"), width=7, height=5, dpi=120,
  ggplot(DT, aes(sd, colour=mbin)) + geom_density() + coord_cartesian(xlim=c(-.1,.1)) +
    theme_minimal() + labs(x="freq_v6 - freq_v3", y="density", colour="minor-allele freq",
    title="v3 vs v6: signed freq-diff density by MAF bin"))
cat("\nplots -> logs/figures/  (sync brings them over)\n")
