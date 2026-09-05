# design_power.R -- what three conventional designs could have detected,
# against this experiment, on one common framework.
#
#   Rscript temp_aging/design_power.R
#
# Designs are Tony's, 4 Sep 2026:
#   A  1000 outbred diploid individuals, one measure each, MAF 0.15, 1e6 tests
#   B  200 inbred lines (haploid chromosomes), 10 measures each, MAF 0.15, 1e6 tests
#   C  750 8-way RILs, 8-way ANOVA on 7 df, 10 measures each, 1e4 tests
#
# All assume h2 = 50% and a locus explaining 2% of individual phenotypic
# variance in the base population. Thresholds are Bonferroni at 0.05 over the
# stated number of tests; this experiment is held to its own LWP thresholds
# instead, which are far stricter than a Bonferroni over 268 tiles would be.
#
# Non-centrality is n * R2 throughout, R2 being the locus contribution as a
# fraction of the variance of whatever unit is measured -- individuals in A,
# line means in B and C.
#
# Base population: Vp = 1 per individual, h2 = 0.5, so Vg = 0.5, Ve = 0.5.
# The locus explains 2% of individual phenotypic variance in that base population.
#
# Inbreeding doubles the locus contribution and the genetic variance: for a
# locus with alleles at frequency p_i and effects a_i, outbred Va = 2*sum p(a-abar)^2
# but among fully inbred lines the variance is sum p(2a - 2abar)^2 / 2 = 2x that.
# Replication cuts only the environmental part: Var(line mean) = Vg' + Ve/r.

Vp <- 1; h2 <- 0.5; Vg <- h2*Vp; Ve <- Vp - Vg; qtl <- 0.02

pw <- function(lam, df, ntest, alpha = 0.05) {
  thr <- qchisq(alpha/ntest, df, lower.tail = FALSE)
  pchisq(thr, df, ncp = lam, lower.tail = FALSE)
}
# n needed for 80% power
n80 <- function(R2, df, ntest, alpha = 0.05) {
  thr <- qchisq(alpha/ntest, df, lower.tail = FALSE)
  lam <- uniroot(function(l) pchisq(thr, df, ncp = l, lower.tail = FALSE) - 0.8,
                 c(1e-6, 5000))$root
  c(lam = lam, n = lam/R2)
}

cat("BASE: Vp = 1, h2 = 50%, locus = 2% of phenotypic variance\n\n")

## A -- 1000 outbred diploids, one measure, 1M tests, 1 df
R2a <- qtl/Vp                       # 0.02 of the measured variance
lamA <- 1000*R2a
cat(sprintf("A  1000 outbred diploids, 1 measure, MAF 0.15, 1e6 tests, 1 df\n"))
cat(sprintf("   R2 = %.3f   lambda = %.1f   threshold = %.1f   POWER = %.3f\n",
            R2a, lamA, qchisq(.05/1e6,1,lower.tail=FALSE), pw(lamA,1,1e6)))
r <- n80(R2a,1,1e6); cat(sprintf("   need n = %.0f individuals for 80% power\n\n", r["n"]))

## B -- 200 inbred lines, 10 measures, 1M tests, 1 df
VgB   <- 2*Vg                       # inbreeding doubles genetic variance
qtlB  <- 2*qtl                      # and the locus contribution
VlineB<- VgB + Ve/10
R2b   <- qtlB/VlineB
lamB  <- 200*R2b
cat(sprintf("B  200 inbred lines, 10 measures each, MAF 0.15, 1e6 tests, 1 df\n"))
cat(sprintf("   locus %.3f of line variance %.3f -> R2 = %.4f\n", qtlB, VlineB, R2b))
cat(sprintf("   lambda = %.1f   threshold = %.1f   POWER = %.3f\n",
            lamB, qchisq(.05/1e6,1,lower.tail=FALSE), pw(lamB,1,1e6)))
r <- n80(R2b,1,1e6); cat(sprintf("   need n = %.0f lines for 80% power\n\n", r["n"]))

## C -- 750 8-way RILs, 10 measures, 10K tests, 7 df
VgC   <- 2*Vg; qtlC <- 2*qtl; VlineC <- VgC + Ve/10
R2c   <- qtlC/VlineC
lamC  <- 750*R2c
cat(sprintf("C  750 8-way RILs, 10 measures each, 1e4 tests, 7 df\n"))
cat(sprintf("   R2 = %.4f   lambda = %.1f   threshold = %.1f   POWER = %.3f\n",
            R2c, lamC, qchisq(.05/1e4,7,lower.tail=FALSE), pw(lamC,7,1e4)))
r <- n80(R2c,7,1e4); cat(sprintf("   need n = %.0f RILs for 80% power\n\n", r["n"]))

## X-QTL as run -- lambda from the fitted relation h2 = 0.00774 (T-7),
## held to the paper's own LWP thresholds rather than a Bonferroni over tiles.
c0 <- 0.00774
lamX <- qtl*100/c0
for (lwp in c(7.5, 15)) {
  thr <- qchisq(10^-lwp, 7, lower.tail = FALSE)
  lam80 <- uniroot(function(l) pchisq(thr,7,ncp=l,lower.tail=FALSE)-0.8, c(1,600))$root
  cat(sprintf("X  this experiment at LWP %4.1f, 7 df\n", lwp))
  cat(sprintf("   lambda at a 2% locus = %.0f   threshold = %.1f   POWER = %.4f\n",
      lamX, thr, pchisq(thr,7,ncp=lamX,lower.tail=FALSE)))
  cat(sprintf("   smallest locus at 80% power = %.2f%%\n\n", c0*lam80))
}
cat("\nSmallest locus detectable at 80% power, each design as specified:\n")
cat(sprintf("  A  %5.2f%% of phenotypic variance\n", 100*n80(1,1,1e6)["lam"]/1000))
cat(sprintf("  B  %5.2f%%\n", 100*n80(1,1,1e6)["lam"]/200*VlineB/2))
cat(sprintf("  C  %5.2f%%\n", 100*n80(1,7,1e4)["lam"]/750*VlineC/2))
cat(sprintf("  X  %5.2f%% (at its LWP 7.5 threshold)\n", c0*54.3))
