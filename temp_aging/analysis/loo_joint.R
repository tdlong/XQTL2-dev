# loo_joint.R -- haplotype-based heritability by the joint model (COJO / GREML).
#
#   Rscript temp_aging/analysis/loo_joint.R
#
# Reads temp_aging/numbers/loo_u.txt.gz, written by loo_genome.R on HPC3.
#
# WHY A JOINT MODEL. Summing h2 over tiles over-counts, because one locus shows
# up at every tile it is in LD with. Clumping (prune + threshold) answers that by
# discarding markers, which throws away real loci and leaves an arbitrary window.
# The field's answer is to keep every marker and divide out the correlation
# between them:
#
#     h2 = u' R^-1 u        instead of        h2 = sum(u^2)
#
# That quadratic form IS SNP-based heritability -- GREML fits exactly this
# through the relationship matrix. Dividing the genome-wide sum by a mean LD
# score is the scalar approximation to it, valid only if R is homogeneous.
#
# WHY IT FACTORISES. The between-tile correlation applies alike to each founder,
# so the full 2128 x 2128 structure is Kronecker, R (x) I_8, and its inverse
# factorises. Each founder's profile along a chromosome is handled on its own
# against a 266 x 266 R, and the results add. No large matrix is ever built.
#
# WHY IT STAYS UNBIASED. u_tr comes from nine cages and u_te from the tenth, so
# their sampling noise is independent and E[u_tr' R^-1 u_te] = u_true' R^-1
# u_true. No bias correction, and -- since tiles are not selected at all here --
# no winner's curse either.
#
# WHERE R COMES FROM. NOT from correlating founder effects pairwise: each entry
# would rest on eight numbers with a sum-to-zero constraint (7 df), and 26% of
# genuinely unlinked pairs exceed |r| = 0.5 at that df. Instead R_ij = rho(d_ij)
# with rho measured as the CROSS-FOLD correlation at lag d -- u_tr(i) against
# u_te(i+d). Train and test noise are independent, so that estimate carries no
# noise contamination at any lag, which the within-sample autocorrelation does
# (the +/-100 kb smoothing correlates neighbouring tiles' errors).

suppressMessages({library(tidyverse); library(Matrix)})

IN    <- "temp_aging/numbers/loo_u.txt.gz"
MAXD  <- 15      # lags used for rho; it reaches the between-arm null by ~10 cM
RIDGE <- 1e-3    # on the correlation scale, so R stays positive definite

if (!file.exists(IN)) stop("missing ", IN,
  "\nRun temp_aging/analysis/loo_genome.sh on HPC3 and pull.", call. = FALSE)
U <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble() %>%
  mutate(trt = paste0(sugar, "_", sex))
cat("tiles:", n_distinct(paste(U$chr, U$tile)), " founders:", n_distinct(U$founder),
    " folds:", n_distinct(U$holdout), " treatments:", n_distinct(U$trt), "\n\n")

# ---- rho(d): signal correlation at lag d, from the cross-fold product --------
# cor(u_tr at i, u_te at i+d) has no noise term at any lag, including d = 0.
rho_hat <- function(z) {
  z <- z %>% arrange(trt, holdout, founder, chr, tile)
  map_dfr(0:MAXD, function(d) {
    w <- z %>% group_by(trt, holdout, founder, chr) %>%
      mutate(b = lead(u_te, d)) %>% ungroup() %>% filter(!is.na(b))
    tibble(d = d, num = mean(w$u_tr * w$b))
  }) %>% mutate(rho = num / num[d == 0])
}

# ---- h2 = sum_f u_f(tr)' R^-1 u_f(te), per arm, summed ----------------------
joint_h2 <- function(z, rho) {
  arms <- z %>% group_split(chr)
  sum(map_dbl(arms, function(a) {
    tl <- sort(unique(a$tile)); n <- length(tl)
    if (n < 2) return(sum(100 * a$k * a$u_tr * a$u_te))
    D <- abs(outer(tl, tl, "-"))
    R <- matrix(approx(rho$d, rho$rho, xout = pmin(D, MAXD), rule = 2)$y, n, n)
    R <- (R + t(R)) / 2; diag(R) <- 1 + RIDGE
    Ri <- tryCatch(solve(R), error = function(e) MASS::ginv(R))
    sum(map_dbl(unique(a$founder), function(f) {
      g <- a %>% filter(founder == f) %>% arrange(tile)
      if (nrow(g) != n) return(NA_real_)
      100 * g$k[1] * as.numeric(t(g$u_tr) %*% Ri %*% g$u_te)
    }), na.rm = TRUE)
  }))
}

jk <- function(x) { n <- length(x); sqrt((n - 1) / n * sum((x - mean(x))^2)) }

cat("rho(d) -- signal correlation between tiles d cM apart (cross-fold, noise-free)\n\n")
R_ALL <- U %>% group_split(trt) %>% map_dfr(~rho_hat(.x) %>% mutate(trt = .x$trt[1]))
print(as.data.frame(R_ALL %>% filter(d <= 10) %>% transmute(d, trt, rho = round(rho, 3)) %>%
  pivot_wider(names_from = trt, values_from = rho)), row.names = FALSE)

cat("\nheritability: naive sum over tiles vs the joint estimate\n")
cat("(the joint value MUST be the smaller -- it removes the over-counting)\n\n")
res <- U %>% group_split(trt, holdout) %>% map_dfr(function(z) {
  rho <- R_ALL %>% filter(trt == z$trt[1])
  tibble(trt = z$trt[1], holdout = z$holdout[1],
         naive = sum(100 * z$k * z$u_tr * z$u_te),
         joint = joint_h2(z, rho))
})
print(as.data.frame(res %>% group_by(trt) %>% summarise(
  naive_sum = round(mean(naive), 1),
  joint_h2  = round(mean(joint), 2),
  jk_se     = round(jk(joint), 2),
  over_count_factor = round(mean(naive) / mean(joint), 2), .groups = "drop")),
  row.names = FALSE)

cat("\nThe over-count factor is the LD score, arrived at without assuming one.\n")
