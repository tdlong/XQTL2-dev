# make_LOO_fig.R -- can nine replicates predict the tenth?
#
#   bash temp_aging/analysis/run_numbers.sh loo_regions   # panel B's numbers
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
# PANEL B. The fraction of marker-based heritability a predictor on the top N
# markers explains in the held-out replicate, with a marker being a single tile
# or a tile and its neighbours within 3 cM. See loo_regions.R for why this
# fraction needs no LD correction while a percent-of-variance version does.

suppressMessages({library(tidyverse); library(patchwork)})

POS   <- "temp_aging/numbers/loo_positions.txt.gz"
REGF  <- "temp_aging/numbers/loo_regions.tsv"  # loo_regions.R
for (f in c(POS, REGF)) if (!file.exists(f)) stop("missing ", f, "\nRun loo_genome.sh on HPC3, pull, then run_numbers.sh loo_regions.", call. = FALSE)
OUT   <- "temp_aging/analysis/LOO_fig.png"
NSHOW <- 30      # positions shown; loo_genome.R clumps at 5 cM
BASE  <- 13
source("temp_aging/analysis/palette.R")   # FOUNDER_COL, TRT_LAB


P <- read.table(POS, header = TRUE, sep = "\t") %>% as_tibble() %>%
  mutate(key = paste0(sugar, "_", sex),
         trt = factor(TRT_LAB[key], levels = TRT_LAB))

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
# Fraction of marker-based heritability explained in the held-out replicate,
# from loo_regions.R. A predictor's marker is either a single tile (W = 0) or
# the tile and its neighbours within 3 cM (W = 3); the dashed line is the
# all-tile predictor. Counted this way the LD over-count cancels (see
# loo_regions.R), so no LD score enters the panel.
REG <- read.table(REGF, header = TRUE, sep = "\t") %>% as_tibble() %>%
  mutate(trt = factor(TRT_LAB[trt], levels = TRT_LAB))
cum <- REG %>% filter(W %in% c(0, 3)) %>% mutate(marker = factor(ifelse(W == 0, "single tile", "tile \u00b1 3 cM"),
                                                             levels = c("single tile", "tile \u00b1 3 cM")))
ALL <- REG %>% filter(is.infinite(W))

pB <- ggplot(cum, aes(N, frac, colour = marker, fill = marker)) +
  geom_hline(data = ALL, aes(yintercept = frac), linetype = "dashed", colour = "grey30") +
  geom_ribbon(aes(ymin = pmax(frac - se, 0), ymax = frac + se), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  facet_wrap(~trt, nrow = 1) +
  scale_x_continuous(breaks = c(1, 10, 20, 30)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_colour_manual(values = c("#1F4E9C", "#C0392B"), name = NULL) +
  scale_fill_manual(values = c("#1F4E9C", "#C0392B"), name = NULL) +
  labs(x = "markers in predictor", y = "fraction of marker-based\nheritability explained") +
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
cat("\npanel B: fraction explained at N = 10 and 30, and by every tile:\n")
print(as.data.frame(REG %>% filter((N %in% c(10, 30) & W %in% c(0, 3)) | is.infinite(W)) %>%
  transmute(trt, W, N, frac = round(frac, 2), se = round(se, 2)) %>% arrange(trt, W, N)), row.names = FALSE)
