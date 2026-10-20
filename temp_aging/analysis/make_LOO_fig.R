# make_LOO_fig.R -- can nine cages predict the tenth?
#
#   Rscript temp_aging/analysis/loo_joint.R      # writes loo_total.txt first
#   Rscript temp_aging/analysis/make_LOO_fig.R
#
# Genome-wide, from loo_genome.R's output. Positions are ranked on the NINE
# training cages and h2 is read in the TENTH, so the held-out cage neither sets
# the weights nor chooses where to look.
#
# PANEL A. The prediction itself. Each point is one founder at one position in
# one held-out cage: x is what the other nine predicted, y is what that cage did.
# The held-out cage's selected proportion is known, so the prediction is
# i_test x e_train with nothing fitted to it -- the dashed line is 1:1, not a
# regression.
#
# PANEL B. The heritability a score on the top N markers explains, against the
# marker-based heritability of the whole genome. Both are percent of phenotypic
# variance, so they are directly comparable and no ratio is needed.
#
# SPACING. Markers stay at loo_genome.R's 5 cM clumping. An earlier version
# re-clumped to 10 cM on the theory that 5 cM double counts; the N x N matrix
# shows it does not -- the naive and matrix-corrected sums agree to ~1% (11.13 vs
# 11.24 at N = 30). The fall seen on re-clumping was losing markers, 30 down to
# 18, not removing redundancy.

suppressMessages({library(tidyverse); library(patchwork)})

POS   <- "temp_aging/numbers/loo_positions.txt.gz"
TOTF  <- "temp_aging/numbers/loo_total.txt"
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
  mutate(trt = factor(TRT_LAB[trt], levels = TRT_LAB))


# ── panel A ──────────────────────────────────────────────────────────────────
A <- P %>% filter(rank <= 10)
lab <- A %>% group_by(trt) %>%
  summarise(r = cor(pred, obs), slope = coef(lm(obs ~ pred))[2], .groups = "drop") %>%
  mutate(txt = sprintf("r = %.2f\nslope = %.2f", r, slope))
L <- max(abs(c(A$pred, A$obs)))

pA <- ggplot(A, aes(pred, obs, colour = founder)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(size = 1.4, alpha = 0.7) +
  geom_text(data = lab, aes(-L * 0.95, L * 0.86, label = txt), inherit.aes = FALSE,
            hjust = 0, size = 3.5, lineheight = 0.95) +
  scale_colour_manual(values = FOUNDER_COL, name = NULL) +
  coord_cartesian(xlim = c(-L, L), ylim = c(-L, L)) +
  facet_wrap(~trt, nrow = 1) +
  labs(x = "predicted frequency change (from nine cages)",
       y = "observed\n(held-out cage)") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3))) +
  theme_bw(BASE) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                         aspect.ratio = 1,
                         strip.background = element_rect(fill = "grey92", colour = NA))

# ── panel B ──────────────────────────────────────────────────────────────────
# What a score on the top N markers explains, against the genome-wide total.
# Both axes are percent of phenotypic variance. The correlation between the N
# markers is divided out with the N x N matrix -- the condition number is 1.4 to
# 3.1 at these spacings, so this is the exact correction, not an approximation.
jk <- function(x) { n <- length(x); sqrt((n - 1) / n * sum((x - mean(x))^2)) }
qform <- function(g, rho) {
  key <- g %>% distinct(rank, chr, cM) %>% arrange(rank); n <- nrow(key)
  R <- ifelse(outer(key$chr, key$chr, "=="),
              rho[pmin(round(abs(outer(key$cM, key$cM, "-"))), length(rho) - 1) + 1], 0)
  diag(R) <- 1; Ri <- solve(R)
  v <- g %>% mutate(v = pred / sqrt(p), w = obs / sqrt(p))
  cc <- sum(g$h2_cross) / sum(v$v * v$w)          # = 100k/i^2, recovered
  q <- function(col) sum(map_dbl(unique(g$founder), function(f) {
    h <- v %>% filter(founder == f) %>% arrange(rank)
    if (nrow(h) != n) return(NA_real_)
    cc * as.numeric(t(h$v) %*% Ri %*% h[[col]]) }), na.rm = TRUE)
  tibble(cross = q("w"), nine = q("v"))
}
RHO <- read.table(RHOF, header = TRUE, sep = "\t") %>% as_tibble()
rl <- split(RHO$rho, RHO$trt)
cum <- map_dfr(1:NSHOW, function(N)
  P %>% filter(rank <= N) %>% group_split(trt, holdout) %>%
    map_dfr(~qform(.x, rl[[as.character(.x$key[1])]]) %>%
              mutate(trt = .x$trt[1], holdout = .x$holdout[1], N = N))) %>%
  group_by(trt, N) %>% summarise(explained = mean(cross)^2 / mean(nine),
    se = jk(cross^2 / nine), .groups = "drop")

pB <- ggplot(cum, aes(N, explained)) +
  geom_hline(data = TOT, aes(yintercept = total), linetype = "dashed",
             colour = "grey30") +
  geom_ribbon(aes(ymin = pmax(explained - se, 0), ymax = explained + se),
              fill = "#1F4E9C", alpha = 0.18) +
  geom_line(colour = "#1F4E9C", linewidth = 0.8) +
  geom_point(colour = "#1F4E9C", size = 1.4) +
  facet_wrap(~trt, nrow = 1) +
  scale_x_continuous(breaks = c(1, 10, 20, 30)) +
  expand_limits(y = 0) +
  labs(x = "markers in score", y = "heritability (%)") +
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
  transmute(trt, N, explained = round(explained, 2)) %>%
  left_join(TOT %>% transmute(trt, total = round(total, 2), se = round(jk_se, 2)),
            by = "trt") %>%
  mutate(frac_of_total = round(explained / total, 2))), row.names = FALSE)
