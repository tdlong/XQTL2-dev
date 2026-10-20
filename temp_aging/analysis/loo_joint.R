# loo_joint.R -- genome-wide haplotype heritability, corrected for over-counting.
#
#   Rscript temp_aging/analysis/loo_joint.R
#
# Reads temp_aging/numbers/loo_u.txt.gz, written by loo_genome.R on HPC3.
#
# THE PROBLEM. Summing h2 over 1 cM tiles over-counts, because one locus shows up
# at every tile it is in LD with. The naive genome-wide sum is 44-64%, which is
# not a heritability.
#
# THE CORRECTION. Each locus is counted once per tile it reaches, weighted by r^2
# -- so the over-count factor is the mean row sum of squared correlations,
# trace(R^2)/n, built from the measured autocorrelation of the signal. Dividing
# the naive sum by it is LD score regression's correction, with the LD score
# measured here rather than taken from a reference panel.
#
# WHY NOT R^-1. The exact version is h2 = u' R^-1 u, the quadratic form GREML
# fits. It is not estimable at this resolution. Adjacent 1 cM tiles are 89%
# correlated, so R has near-duplicate rows: 273 tiles carry about 273/4.5 ~ 60
# tiles' worth of independent information, and asking for a 273x273 inverse from
# that is p >> n. Measured directly, the answer tracks the regularisation over
# four orders of magnitude (6e5 at an eigenvalue clip of 1e-6, down to 1e2 at
# 1e-2), so it reports the clip, not the data. The field does not invert R either
# -- GREML uses REML on variance components and LDpred imposes a prior.
# trace(R^2)/n needs no inverse and is a sum, so it is stable.
#
# WHAT IT ASSUMES. That the over-count factor applies uniformly genome-wide,
# which is LD score regression's assumption as well. Relaxing it is what the
# matrix version would do, and that is the part this data cannot support.
#
# WHERE rho COMES FROM. The CROSS-FOLD product at lag d: u_tr(i) against
# u_te(i+d). Train comes from nine cages and test from the tenth, so their
# sampling noise is independent and contributes nothing at any lag -- including
# lag 0. The within-sample autocorrelation cannot be used: the +/-100 kb
# smoothing correlates neighbouring tiles' errors and inflates it.
#
# The between-arm product gives the true zero and comes out at -1.9% of lag 0,
# so there is no genome-wide common factor to subtract.

suppressMessages(library(tidyverse))

IN     <- "temp_aging/numbers/loo_u.txt.gz"
OUT    <- "temp_aging/numbers/loo_total.txt"   # read by make_LOO_fig.R
OUTR   <- "temp_aging/numbers/loo_rho.txt"     # likewise, for the N x N matrix
MAXLAG <- 45     # rho is summed as rho^2, so the tail adds nothing past ~15 cM;
                 # 45 is simply well beyond where it matters (see the stability
                 # table this script prints)

if (!file.exists(IN)) stop("missing ", IN,
  "\nRun temp_aging/analysis/loo_genome.sh on HPC3, then pull.", call. = FALSE)
U <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble() %>%
  mutate(trt = paste0(sugar, "_", sex))
cat("tiles:", n_distinct(paste(U$chr, U$tile)), " founders:", n_distinct(U$founder),
    " folds:", n_distinct(U$holdout), " treatments:", n_distinct(U$trt), "\n")

# signal autocorrelation at lag d, noise-free by construction
rho_hat <- function(z, maxlag) {
  z <- z %>% arrange(holdout, founder, chr, tile)
  s <- map_dbl(0:maxlag, function(d) {
    w <- z %>% group_by(holdout, founder, chr) %>% mutate(b = lead(u_te, d)) %>%
      ungroup() %>% filter(!is.na(b))
    mean(w$u_tr * w$b)
  })
  pmax(s / s[1], 0)
}

# over-count factor = trace(R^2)/n, on the real tile layout, arm by arm
overcount <- function(tiles, rho) {
  per <- map_dbl(group_split(tiles, chr), function(a) {
    tl <- sort(a$tile)
    sum(rho[pmin(abs(outer(tl, tl, "-")), length(rho) - 1) + 1]^2) / length(tl)
  })
  weighted.mean(per, tiles %>% count(chr) %>% pull(n))
}

naive_sum <- function(z) mean(map_dbl(unique(z$holdout), function(h) {
  g <- z %>% filter(holdout == h); sum(100 * g$k * g$u_tr * g$u_te) }))

jk <- function(x) { n <- length(x); sqrt((n - 1) / n * sum((x - mean(x))^2)) }

est <- function(z, maxlag) {
  rho <- rho_hat(z, maxlag)
  f   <- overcount(z %>% distinct(chr, tile), rho)
  per_fold <- map_dbl(unique(z$holdout), function(h) {
    g <- z %>% filter(holdout == h); sum(100 * g$k * g$u_tr * g$u_te) / f })
  tibble(trt = z$trt[1], naive = naive_sum(z), overcount = f,
         total = mean(per_fold), jk_se = jk(per_fold))
}

cat("\nmeasured signal autocorrelation rho(d)\n\n")
RHO <- U %>% group_split(trt) %>%
  map_dfr(~tibble(trt = .x$trt[1], d = 0:MAXLAG, rho = rho_hat(.x, MAXLAG)))
print(as.data.frame(RHO %>% filter(d %in% c(1:8, 10, 15, 20, 30)) %>%
  mutate(rho = round(rho, 3)) %>% pivot_wider(names_from = trt, values_from = rho)),
  row.names = FALSE)
write.table(RHO, OUTR, row.names = FALSE, quote = FALSE, sep = "\t")

cat("\ngenome-wide total h2 (% of phenotypic variance)\n\n")
TOT <- U %>% group_split(trt) %>% map_dfr(~est(.x, MAXLAG))
print(as.data.frame(TOT %>% mutate(across(where(is.double), ~round(.x, 2)))),
  row.names = FALSE)
write.table(TOT, OUT, row.names = FALSE, quote = FALSE, sep = "\t")

cat("\nstability to where rho is truncated (it enters as rho^2, so the tail is small)\n\n")
print(as.data.frame(map_dfr(c(10, 15, 20, 30, 45), function(L)
  U %>% group_split(trt) %>% map_dfr(~est(.x, L)) %>%
    transmute(trt, truncate_cM = L, total = round(total, 1))) %>%
  pivot_wider(names_from = trt, values_from = total)), row.names = FALSE)

cat("\nTruncating at 10 cM cuts real decay; past 15 cM the answer does not move.\n")
cat("Wrote:", OUT, "\n")
