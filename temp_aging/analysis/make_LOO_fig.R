# make_LOO_fig.R -- can nine replicates predict the tenth?
#
#   bash temp_aging/analysis/run_numbers.sh yardstick   # marker-based h2, the dashed line
#   Rscript temp_aging/analysis/loo_joint.R             # rho(d), for the N x N matrix
#   Rscript temp_aging/analysis/make_LOO_fig.R
#
# Genome-wide, from loo_genome.R's output. Positions are ranked on the NINE
# replicates and h2 is read in the held-out TENTH, so the held-out replicate
# neither sets the weights nor chooses where to look.
#
# PANEL A. The prediction itself. Each point is one founder at one position in
# one held-out cage: x is what the other nine predicted, y is what that cage did.
# The held-out cage's selected proportion is known, so the prediction is
# i_test x e_train with nothing fitted to it -- the dashed line is 1:1, not a
# regression.
#
# PANEL B. The heritability a predictor on the top N markers explains, against
# the marker-based heritability of the whole genome. Both are percent of
# phenotypic variance. Uncertainty is shown as the pieces themselves: the ten
# held-out replicates as light curves, the ten leave-one-replicate-out
# heritabilities as light dashes.
#
# SPACING. Markers stay at loo_genome.R's 5 cM clumping; the N x N correlation
# matrix divides out what correlation remains (up to 5% of the sum).

suppressMessages({library(tidyverse); library(patchwork)})

POS   <- "temp_aging/numbers/loo_positions.txt.gz"
TOTF  <- "temp_aging/numbers/yardstick.tsv"   # marker-based h2, yardstick.R
PKF   <- "temp_aging/numbers/yardstick_peaks.tsv" # the per-peak divisors behind it
PKF   <- "temp_aging/numbers/yardstick_peaks.tsv" # the per-peak divisors behind it
RHOF  <- "temp_aging/numbers/loo_rho.txt"
OUT   <- "temp_aging/analysis/LOO_fig.png"
NSHOW <- 30      # positions shown; loo_genome.R clumps at 5 cM
BASE  <- 13
source("temp_aging/analysis/palette.R")   # FOUNDER_COL, TRT_LAB

for (f in c(POS, TOTF, RHOF)) if (!file.exists(f)) stop("missing ", f,
  "\nRun loo_genome.sh on HPC3, pull, then loo_joint.R.", call. = FALSE)

P <- read.table(POS, header = TRUE, sep = "\t") %>% as_tibble() %>%
  mutate(key = paste0(sugar, "_", sex),
         trt = factor(TRT_LAB[key], levels = TRT_LAB))
TOT <- read.table(TOTF, header = TRUE, sep = "\t") %>% as_tibble() %>%
  rename(total = h2, jk_se = se) %>% mutate(trt = factor(TRT_LAB[trt], levels = TRT_LAB))
PK  <- read.table(PKF, header = TRUE, sep = "\t") %>% as_tibble()
PK  <- read.table(PKF, header = TRUE, sep = "\t") %>% as_tibble()


# ── panel A ──────────────────────────────────────────────────────────────────
A <- P %>% filter(rank <= 10)
lab <- A %>% group_by(trt) %>%
  summarise(r = cor(pred, obs), slope = coef(lm(obs ~ pred))[2], .groups = "drop") %>%
  mutate(txt = sprintf("R\u00b2 = %.2f", r^2))
L <- max(abs(c(A$pred, A$obs)))

pA <- ggplot(A, aes(pred, obs, colour = founder)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(size = 1.4, alpha = 0.7) +
  geom_text(data = lab, aes(-L * 0.95, L * 0.86, label = txt), inherit.aes = FALSE,
            hjust = 0, size = 3.5, lineheight = 0.95) +
  scale_colour_manual(values = FOUNDER_COL, name = NULL) +
  coord_cartesian(xlim = c(-L, L), ylim = c(-L, L)) +
  facet_wrap(~trt, nrow = 1) +
  labs(x = "predicted frequency change (from nine replicates)",
       y = "observed\n(held-out replicate)") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3))) +
  theme_bw(BASE) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                         aspect.ratio = 1,
                         strip.background = element_rect(fill = "grey92", colour = NA))

# ── panel B ──────────────────────────────────────────────────────────────────
# Heritability a predictor on the N best tiles explains in the held-out
# replicate, as a fraction of phenotypic variance.  Per fold:
#   cross = sum over the N tiles of  100 k sum_f  e_bar_f e_f(r) / p_f   (= sum of h2_r)
#   nine  = sum over the N tiles of  100 k sum_f  e_bar_f e_bar_f / p_f  (the predictor's variance)
#   explained = mean(cross)^2 / mean(nine)   over the ten folds
# which is the regression of the held-out replicate on the predictor, cov^2/var,
# in h2 units.  Plain sums over tiles: the tiles are >= 5 cM apart.  The N x N
# correlation matrix built from rho(d) is computed alongside and printed as a
# check; it moves the sums by at most 5% and is not what is plotted.
jk <- function(x) { n <- length(x); sqrt((n - 1) / n * sum((x - mean(x))^2)) }
sums <- function(g, rho) {
  key <- g %>% distinct(rank, chr, cM) %>% arrange(rank); n <- nrow(key)
  R <- ifelse(outer(key$chr, key$chr, "=="),
              rho[pmin(round(abs(outer(key$cM, key$cM, "-"))), length(rho) - 1) + 1], 0)
  diag(R) <- 1; Ri <- solve(R)
  v <- g %>% mutate(v = pred / sqrt(p), w = obs / sqrt(p))
  cc <- sum(g$h2_cross) / sum(v$v * v$w)          # = 100k/i^2, recovered
  q <- function(col, M) sum(map_dbl(unique(g$founder), function(f) {
    h <- v %>% filter(founder == f) %>% arrange(rank)
    if (nrow(h) != n) return(NA_real_)
    cc * as.numeric(t(h$v) %*% M %*% h[[col]]) }), na.rm = TRUE)
  tibble(cross = q("w", diag(n)), nine = q("v", diag(n)),
         cross_C = q("w", Ri), nine_C = q("v", Ri))
}
RHO <- read.table(RHOF, header = TRUE, sep = "\t") %>% as_tibble()
rl <- split(RHO$rho, RHO$trt)
folds <- map_dfr(1:NSHOW, function(N)
  P %>% filter(rank <= N) %>% group_split(trt, holdout) %>%
    map_dfr(~sums(.x, rl[[as.character(.x$key[1])]]) %>%
              mutate(trt = .x$trt[1], holdout = .x$holdout[1], N = N))) %>%
  mutate(fold = cross_C^2 / nine_C)                # each held-out replicate on its own
cum <- folds %>% group_by(trt, N) %>%
  summarise(explained = mean(cross_C)^2 / mean(nine_C),
            plain_sum = mean(cross)^2 / mean(nine), .groups = "drop")

# The divisor's spread: each Table S1 peak's footprint gives its own tiles per
# locus, so its own genome-wide total = tile sum / footprint. Light grey band:
# the 25-75% quantiles of those totals. Dashed: the pooled estimate.
BAND <- TOT %>% group_split(trt) %>% map_dfr(function(z) {
  v <- z$sum_tiles / PK$tiles_per_locus
  tibble(trt = z$trt, lo = quantile(v, .25), hi = quantile(v, .75)) })

# light blue: the ten held-out replicates one at a time; dark: all ten combined
pB <- ggplot(cum, aes(N, explained)) +
  geom_rect(data = BAND, aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi),
            fill = "grey88", inherit.aes = FALSE) +
  geom_hline(data = TOT, aes(yintercept = total), linetype = "dashed", colour = "grey20") +
  geom_line(data = folds, aes(N, fold, group = holdout), colour = "#1F4E9C", alpha = 0.25, linewidth = 0.4) +
  geom_line(colour = "#1F4E9C", linewidth = 0.9) +
  geom_point(colour = "#1F4E9C", size = 1.4) +
  facet_wrap(~trt, nrow = 1) +
  scale_x_continuous(breaks = c(1, 10, 20, 30)) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  labs(x = "markers in predictor", y = "heritability (%)") +
  theme_bw(BASE) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                         aspect.ratio = 1,
                         strip.background = element_rect(fill = "grey92", colour = NA))

png(OUT, width = 11, height = 7.8, units = "in", res = 300)
print(pA / pB + plot_layout(guides = "collect") +
      plot_annotation(tag_levels = "A") &
      theme(legend.position = "bottom", legend.box = "vertical",
            legend.margin = margin(0, 0, 0, 0)))
invisible(dev.off())
cat("Wrote:", OUT, "\n\n")

cat("panel A, top 10 positions:\n")
print(as.data.frame(lab %>% transmute(trt, r = round(r, 3),
  slope = round(slope, 3))), row.names = FALSE)
cat("\npanel B endpoints, and the genome-wide total:\n")
print(as.data.frame(cum %>% filter(N %in% c(10, 30)) %>%
  transmute(trt, N, explained = round(explained, 2), plain_sum = round(plain_sum, 2)) %>%
  left_join(TOT %>% transmute(trt, total = round(total, 2)),
            by = "trt") %>%
  mutate(frac_of_total = round(explained / total, 2))), row.names = FALSE)
