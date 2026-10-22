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
CURF  <- "temp_aging/numbers/loo_curves.tsv"   # both curves, loo_curves.R
PKF   <- "temp_aging/numbers/yardstick_peaks.tsv" # the per-peak divisors behind it
CURF  <- "temp_aging/numbers/loo_curves.tsv"   # both curves, loo_curves.R
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
# Two curves against the marker-based heritability of the trait (dashed), all in
# percent of phenotypic variance (loo_curves.R):
#   CEILING, how much heritability the top N regions hold -- ranked in five
#     replicates, effects from four others, scored on the held-out one, so no
#     winner's curse. Conservative: the predictor ranks on nine.
#   PREDICTOR, what an N-marker predictor built on the nine explains in the
#     held-out replicate, with the ten replicates drawn individually in light blue.
# The gap up to the dashed line is heritability outside the N regions; the gap
# between the curves is what is lost to estimating the effects.
CUR <- read.table(CURF, header = TRUE, sep = "\t") %>% as_tibble() %>%
  mutate(trt = factor(TRT_LAB[trt], levels = TRT_LAB))
folds <- map_dfr(1:max(CUR$N), function(N) P %>% filter(rank <= N) %>% group_by(trt, holdout) %>%
  summarise(fold = sum(h2_cross)^2 / sum(h2_nine), .groups = "drop") %>% mutate(N = N))
BAND <- TOT %>% group_split(trt) %>% map_dfr(function(z) {
  v <- z$sum_tiles / PK$tiles_per_locus
  tibble(trt = z$trt, lo = quantile(v, .25), hi = quantile(v, .75)) })
LONG <- CUR %>% select(trt, N, ceiling, predictor) %>%
  pivot_longer(c(ceiling, predictor), names_to = "what", values_to = "h2") %>%
  mutate(what = factor(ifelse(what == "ceiling", "heritability at the N regions",
                              "explained by the N-marker predictor"),
                       levels = c("heritability at the N regions", "explained by the N-marker predictor")))

pB <- ggplot(LONG, aes(N, h2, colour = what)) +
  geom_rect(data = BAND, aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi),
            fill = "grey88", inherit.aes = FALSE) +
  geom_hline(data = TOT, aes(yintercept = total), linetype = "dashed", colour = "grey20") +
  geom_line(data = folds, aes(N, fold, group = holdout), colour = "#1F4E9C",
            alpha = 0.3, linewidth = 0.4, inherit.aes = FALSE) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("heritability at the N regions" = "#C0392B",
                                 "explained by the N-marker predictor" = "#1F4E9C"), name = NULL) +
  facet_wrap(~trt, nrow = 1) +
  scale_x_continuous(breaks = c(1, 10, 20, 30)) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  labs(x = "markers in predictor", y = "heritability (%)") +
  guides(colour = guide_legend(nrow = 1)) +
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
cat("\npanel B at N = 10, 30 and the end:\n")
print(as.data.frame(CUR %>% filter(N %in% c(10, 30, max(N))) %>%
  transmute(trt, N, `regions hold` = round(ceiling, 1), `predictor explains` = round(predictor, 1),
            `trait h2` = round(trait, 1))), row.names = FALSE)
