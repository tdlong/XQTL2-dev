# design_power.R -- what three conventional designs could have detected.
#
#   Rscript temp_aging/design_power.R
#
# Designs are Tony's, 4 Sep 2026. Autosomes only.
#   A  1000 outbred diploid individuals, one measure each, MAF 0.15, 1e6 tests
#   B  200 inbred lines (haploid chromosomes), 10 measures each, MAF 0.15, 1e6 tests
#   C  750 8-way RILs, 8-way ANOVA on 7 df, 10 measures each, 1e4 tests
#
# All assume h2 = 50% and a locus explaining 2% of individual phenotypic
# variance in the base population. Thresholds are Bonferroni at 0.05 over the
# stated number of tests.
#
# lambda = n * R2 / (1 - R2), R2 being the locus contribution as a fraction of
# the variance of whatever unit is measured -- individuals in A, line means in
# B and C. The script derives every quantity rather than asserting it, then
# checks the analytic power against simulation.
#
# THE POINT THAT IS EASY TO GET WRONG. Inbreeding doubles the locus variance,
# but it doubles the unlinked background too, so the two cancel in the ratio.
# Locus R2 goes 2.00% -> 2.67% from homozygosity alone; it only reaches 3.81%
# because ten measures per line cut Ve tenfold. Nearly all of the inbred-line
# advantage is the replication, not the homozygosity.

set.seed(1)
p <- 0.15; q <- 1-p          # p = frequency of the minor / 'a' allele
Vp <- 1; h2 <- 0.5
Vg <- h2*Vp; Ve <- Vp - Vg   # 0.5 and 0.5, per individual, outbred base
QTL <- 0.02                  # locus explains 2% of individual phenotypic variance
r <- 10                      # measures per line in designs B and C

cat("=== base population, outbred diploid ===\n")
cat(sprintf("  genotype freqs  aa %.4f   Aa %.4f   AA %.4f\n", p^2, 2*p*q, q^2))
alpha <- sqrt(QTL/(2*p*q))   # additive allele effect; Va = 2pq alpha^2
cat(sprintf("  Va(locus) = 2pq*alpha^2 = %.4f  ->  alpha = %.4f\n", 2*p*q*alpha^2, alpha))
cat(sprintf("  Vg total %.2f = locus %.2f + background %.2f ;  Ve %.2f ;  Vp %.2f\n\n",
            Vg, QTL, Vg-QTL, Ve, Vp))

cat("=== inbred lines: only aa and AA, at 15% and 85% ===\n")
vloc_line <- p*q*(2*alpha)^2
cat(sprintf("  locus var among lines = pq(2alpha)^2 = %.4f   (= %.1fx the outbred %.4f)\n",
            vloc_line, vloc_line/QTL, QTL))
vbg_line <- 2*(Vg-QTL)
cat(sprintf("  background Vg among lines = 2 x %.2f = %.2f\n", Vg-QTL, vbg_line))
cat(sprintf("  total genetic among lines = %.2f (= 2 x Vg, as it must be)\n", vloc_line+vbg_line))
cat(sprintf("  %d measures per line -> Ve on a line mean = %.2f/%d = %.3f\n", r, Ve, r, Ve/r))
Vline <- vloc_line + vbg_line + Ve/r
R2_line <- vloc_line/Vline
cat(sprintf("  Var(line mean) = %.3f ;  locus R2 = %.4f\n", Vline, R2_line))
cat(sprintf("  unreplicated for comparison: R2 = %.4f/%.3f = %.4f\n\n",
            vloc_line, vloc_line+vbg_line+Ve, vloc_line/(vloc_line+vbg_line+Ve)))

pw <- function(R2, n, df, ntest) {
  lam <- n*R2/(1-R2)
  thr <- qchisq(0.05/ntest, df, lower.tail=FALSE)
  c(lambda=lam, thr=thr, power=pchisq(thr, df, ncp=lam, lower.tail=FALSE))
}
n80 <- function(R2, df, ntest) {
  thr <- qchisq(0.05/ntest, df, lower.tail=FALSE)
  lam <- uniroot(function(l) pchisq(thr,df,ncp=l,lower.tail=FALSE)-0.8, c(1e-6,5000))$root
  c(lam=lam, n=lam/(R2/(1-R2)))
}

D <- list(A=list(n=1000, R2=QTL/Vp,  df=1, nt=1e6, lab="1000 outbred diploids, 1 measure"),
          B=list(n=200,  R2=R2_line, df=1, nt=1e6, lab="200 inbred lines, 10 measures"),
          C=list(n=750,  R2=R2_line, df=7, nt=1e4, lab="750 8-way RILs, 10 measures"))
cat("=== analytic ===\n")
for (k in names(D)) { d <- D[[k]]; z <- pw(d$R2,d$n,d$df,d$nt); m <- n80(d$R2,d$df,d$nt)
  cat(sprintf("  %s  %-34s n=%4d df=%d  R2=%.4f  lambda=%6.2f  thr=%.1f  POWER=%.4f\n",
      k, d$lab, d$n, d$df, d$R2, z["lambda"], z["thr"], z["power"]))
  cat(sprintf("       need n = %.0f for 80%% power; smallest locus at n as given = %.2f%% of Vp\n",
      m["n"], 100*m["lam"]/d$n*(if (k=="A") 1 else Vline/2))) }

cat("\n=== simulation check, 2000 reps each ===\n")
NS <- 2000
thr1 <- qchisq(0.05/1e6, 1, lower.tail=FALSE)
hitA <- sum(replicate(NS, { g <- rbinom(1000,2,p)
  y <- g*alpha + rnorm(1000,0,sqrt(Vg-QTL)) + rnorm(1000,0,sqrt(Ve))
  summary(lm(y~g))$fstatistic[1] > thr1 }))
hitB <- sum(replicate(NS, { gl <- rbinom(200,1,p)
  y <- gl*(2*alpha) + rnorm(200,0,sqrt(vbg_line)) + rnorm(200,0,sqrt(Ve/r))
  summary(lm(y~gl))$fstatistic[1] > thr1 }))
fe <- scale(rnorm(8))[,1]; fe <- fe*sqrt(vloc_line/mean(fe^2))
thr7 <- qchisq(0.05/1e4, 7, lower.tail=FALSE)
hitC <- sum(replicate(NS, { f <- sample(1:8,750,replace=TRUE)
  y <- fe[f] + rnorm(750,0,sqrt(vbg_line)) + rnorm(750,0,sqrt(Ve/r))
  7*anova(lm(y~factor(f)))[1,4] > thr7 }))
for (k in names(D)) { d <- D[[k]]
  cat(sprintf("  %s  analytic %.4f   simulated %.4f\n", k,
      pw(d$R2,d$n,d$df,d$nt)["power"], c(hitA,hitB,hitC)[match(k,names(D))]/NS)) }
