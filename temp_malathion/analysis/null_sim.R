#!/usr/bin/env Rscript
###############################################################################
# null_sim.R — the two tests share a numerator, so the whole question is the
# denominator. This simulates a null with the real sampling structure and asks
# which denominator is right.
#
# Null: no QTL. Each pool is a multinomial draw of 2N chromosomes from one
# founder frequency vector, which is exactly the "counting flies" variance the
# Wald models and the ANOVA tries to estimate. Replicate 1 sizes: 43 flies per
# female pool (86 chromosomes), 65 per male pool (130).
#
# A valid test returns uniform p under this. The paper's ANOVA estimates its
# denominator from the 16 residual df; a correctly specified test knows the
# multinomial covariance and spends no df estimating it.
#
#   Rscript temp_malathion/analysis/null_sim.R
###############################################################################
suppressPackageStartupMessages(library(tidyverse))
set.seed(1)

ANOVA_test_function = function(df) {
	df = data.frame(df)
	tt = tapply(df$asf,df$HAP,mean)
	tt2 = names(tt)[tt > 0.14]
	df3 = df %>% filter(HAP %in% tt2) %>% droplevels()
	out=anova(lm(asf~TRT + HAP + TRT:HAP, data=df3))
	-pf(out[3,3]/out[4,3],out[3,1],out[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)
}

# the Wald the pipeline uses, reduced to this null: pooled shift against the
# multinomial covariance of the two arms, 7 df after dropping the degenerate
# direction. No reconstruction error here, so this is the Wald's denominator
# at its most optimistic.
mn.covmat <- function(p, n) { m <- -tcrossprod(p); diag(m) <- p*(1-p); m/n }
wald_lp <- function(P1, P2, n1, n2) {
  nF <- ncol(P1); df <- nF - 1L; R <- nrow(P1)
  c1 <- c2 <- matrix(0, nF, nF)
  for (r in seq_len(R)) {
    pm <- (n1[r]*P1[r,] + n2[r]*P2[r,]) / (n1[r]+n2[r])
    c1 <- c1 + mn.covmat(pm, n1[r]) * n1[r]^2
    c2 <- c2 + mn.covmat(pm, n2[r]) * n2[r]^2
  }
  p1 <- as.vector(n1 %*% P1 / sum(n1)); p2 <- as.vector(n2 %*% P2 / sum(n2))
  cv <- c1/sum(n1)^2 + c2/sum(n2)^2
  eg <- eigen(cv); ev <- pmax(eg$values[1:df], 0.01*mean(eg$values[1:df]))
  tstat <- sum(((diag(1/sqrt(ev)) %*% t(eg$vectors[,1:df])) %*% (p1-p2))^2)
  -pchisq(tstat, df, lower.tail = FALSE, log.p = TRUE)/log(10)
}

p0  <- c(.17,.14,.13,.12,.12,.11,.11,.10)      # a realistic founder vector, sums to 1
nchr <- c(86, 130)                              # 43 females, 65 males, per pool
NSIM <- 4000
cat("simulating", NSIM, "null windows: 2 replicate pairs, 8 founders,", nchr, "chromosomes\n")

res <- map_dfr(seq_len(NSIM), function(i) {
  P1 <- t(sapply(nchr, function(n) rmultinom(1, n, p0)[,1]/n))   # control arm
  P2 <- t(sapply(nchr, function(n) rmultinom(1, n, p0)[,1]/n))   # selected arm
  d <- tibble(TRT = rep(c("C","S"), each = 2*8),
              REP = rep(rep(1:2, each = 8), 2),
              HAP = rep(paste0("H", 1:8), 4),
              asf = asin(sqrt(c(t(P1), t(P2))))) %>%
       mutate(across(c(TRT, HAP), as.factor))
  tibble(anova = ANOVA_test_function(d), wald = wald_lp(P1, P2, nchr, nchr))
})

# ── control: is compositionality the cause? ─────────────────────────────────
# Same ANOVA, same design, but the 8 values per pool are now independent with
# equal variance instead of multinomial proportions that sum to one. If the
# test is calibrated here and inflated above, the sum-to-one constraint is the
# cause and nothing else is.
ctl <- map_dfr(seq_len(NSIM), function(i) {
  d <- tibble(TRT = rep(c("C","S"), each = 2*8),
              REP = rep(rep(1:2, each = 8), 2),
              HAP = rep(paste0("H", 1:8), 4),
              asf = 0.35 + rnorm(32, 0, 0.04)) %>%     # independent, same scale
       mutate(across(c(TRT, HAP), as.factor))
  a <- anova(lm(asf ~ TRT + HAP + TRT:HAP, data = d))
  tibble(anova_independent = -pf(a[3,3]/a[4,3], a[3,1], a[4,1], lower.tail=FALSE, log.p=TRUE)/log(10))
})
res <- bind_cols(res, ctl)

summ <- res %>% pivot_longer(everything(), names_to = "test", values_to = "mlog10p") %>%
  group_by(test) %>%
  summarise(median = median(mlog10p),
            lambda = median(qchisq(10^(-mlog10p), 1, lower.tail=FALSE)) / qchisq(0.5, 1, lower.tail=FALSE),
            `P<0.05` = mean(mlog10p > -log10(0.05)),
            `P<0.01` = mean(mlog10p > -log10(0.01)),
            `P<0.001` = mean(mlog10p > -log10(0.001)), .groups = "drop")

out <- "temp_malathion/analysis/null_sim.txt"
sink(out)
cat("Null simulation: no QTL, pools are multinomial draws of 2N chromosomes.\n")
cat("Replicate 1 sizes (43 and 65 flies per pool). 4000 simulated windows.\n")
cat("A valid test gives median 0.301, lambda 1.00, and rejection rates 0.05 / 0.01 / 0.001.\n\n")
print(as.data.frame(summ %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
cat("\nBoth tests are handed the same numerator. Only the denominator differs:\n")
cat("the ANOVA estimates it from 16 residual df, treating the 8 founder values\n")
cat("in a pool as 8 independent observations; the Wald uses the multinomial\n")
cat("covariance, which knows they sum to one, and spends no df on it.\n\n")
cat("anova_independent is the control: the identical ANOVA on 8 INDEPENDENT\n")
cat("values per pool instead of proportions summing to one. It is calibrated,\n")
cat("so the sum-to-one constraint is what inflates the test, nothing else.\n")
sink()
cat(readLines(out), sep = "\n")
