# make_figure_rr.R — recombination rate, on the same axis as Figures 1 and 3.
#
#   Rscript temp_aging/make_figure_rr.R
#
# Why: the Wald scans show a broad elevated region at the base of 3L and the base
# of 3R. Two candidate explanations -- the 75 kb window and 100 kb smoothing
# spreading one signal, or genuinely low recombination over a wide stretch, which
# would put many megabases in linkage and make a real signal that broad.
# This draws the second so it can be compared against the first, with the Wald
# scan averaged over the four treatments on the same axis. Two y scales: cM/Mb on
# the left, mean -log10 P on the right, the second rescaled onto the first purely
# so both fit one panel -- the right-hand axis carries its own numbers and the two
# have no common unit.
#
# Rate is d(cM)/d(Mb) from pipeline/helpfiles/flymap.r6.txt, the same map the
# scans use for cM. That map has ~119 markers per arm, median spacing 271 kb, so
# it cannot resolve anything finer than a few hundred kb whatever smoothing is
# applied. BANDWIDTH is the ksmooth width used before differentiating; the
# pipeline's own add_genetic() uses 3 Mb, and that is the default here so the
# rate shown is the rate the pipeline's cM values imply. Lower it to see finer
# structure, at the cost of noise the map cannot really support.
#
# Grey bands and dotted dividers are Figure 1's: dm6 euchromatin boundaries
# (Huynh et al. 2023 PLoS Genet 19:e1010439, Table S2).

suppressMessages({library(tidyverse); library(patchwork)})

FLYMAP    <- "pipeline/helpfiles/flymap.r6.txt"
SCAN      <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"   # for the axis extents
OUT       <- "temp_aging/FigureRR_plot.png"
BANDWIDTH <- 3e6
STEP      <- 1e4
W_IN <- 7.5; H_IN <- 3.0; DPI <- 300

CHRS   <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
CHRLAB <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")
HET <- tribble(~chr, ~eu_start, ~eu_end,
  "chrX", 2.5, 21.2, "chr2L", 0.5, 22.9, "chr2R", 1.3, 25.1,
  "chr3L", 0.7, 24.0, "chr3R", 4.5, 32.0)

# ── axis extents taken from the scan, so this overlays Figures 1 and 3 ───────
sc <- read.table(SCAN, header = TRUE, sep = "\t") %>% as_tibble() %>%
  filter(chr %in% CHRS)

# mean Wald over the four treatments, per window
wald <- sc %>% group_by(chr, pos) %>%
  summarise(w = mean(Wald_log10p, na.rm = TRUE), .groups = "drop")
lens <- sc %>% group_by(chr) %>% summarise(len = max(pos), .groups = "drop") %>%
  mutate(chr = factor(chr, levels = CHRS)) %>% arrange(chr) %>%
  mutate(offset = lag(cumsum(len), default = 0))

# ── rate = derivative of the smoothed map ────────────────────────────────────
fm <- read.table(FLYMAP, header = FALSE)
colnames(fm) <- c("chr", "pos", "cM")

rr <- lens %>% pmap_dfr(function(chr, len, offset) {
  x <- fm %>% filter(chr == !!chr) %>% arrange(pos)
  o <- ksmooth(x$pos, x$cM, kernel = "normal", bandwidth = BANDWIDTH)
  f <- splinefun(o$x, o$y)
  at <- seq(min(x$pos), min(max(x$pos), len), by = STEP)
  tibble(chr = chr, pos = at,
         # splinefun's first derivative is cM per bp; x1e6 for cM/Mb
         rate = f(at, deriv = 1) * 1e6)
}) %>% mutate(chr = factor(chr, levels = CHRS)) %>%
  left_join(lens %>% select(chr, offset), by = "chr") %>%
  mutate(gx = (pos + offset) / 1e6,
         rate = pmax(rate, 0))          # the spline can dip below 0 between markers

het_bands <- HET %>% mutate(chr = factor(chr, levels = CHRS)) %>%
  left_join(lens, by = "chr") %>%
  { bind_rows(
      transmute(., xmin = offset/1e6, xmax = (offset + eu_start*1e6)/1e6),
      transmute(., xmin = (offset + eu_end*1e6)/1e6, xmax = (offset + len)/1e6)) } %>%
  filter(xmax > xmin)
chr_breaks <- (lens$offset + lens$len/2) / 1e6
chr_edges  <- (lens$offset + lens$len)[-nrow(lens)] / 1e6
xmax_all   <- sum(lens$len) / 1e6

wald <- wald %>% mutate(chr = factor(chr, levels = CHRS)) %>%
  left_join(lens %>% select(chr, offset), by = "chr") %>%
  mutate(gx = (pos + offset) / 1e6)

# right axis rescaled onto the left so both curves share one panel
RMAX <- max(rr$rate); WMAX <- max(wald$w)
wald <- wald %>% mutate(y = w * RMAX / WMAX)

cat("recombination rate, cM/Mb  (ksmooth bandwidth", BANDWIDTH/1e6, "Mb)\n\n")
rr %>% group_by(chr) %>%
  summarise(median = round(median(rate), 2), max = round(max(rate), 2),
            `Mb below 0.5` = round(sum(rate < 0.5) * STEP / 1e6, 1),
            `Mb below 1.0` = round(sum(rate < 1.0) * STEP / 1e6, 1),
            .groups = "drop") %>% as.data.frame() %>% print(row.names = FALSE)
cat("\n")

p <- ggplot(rr) +
  geom_rect(data = het_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey88", colour = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = chr_edges, linetype = "dotted",
             colour = "grey35", linewidth = 0.3) +
  geom_line(data = wald, aes(gx, y, group = chr),
            linewidth = 0.3, colour = "#D62728") +
  geom_line(aes(gx, rate, group = chr), linewidth = 0.5, colour = "grey15") +
  scale_x_continuous(expand = expansion(0), breaks = chr_breaks,
                     labels = CHRLAB[CHRS]) +
  scale_y_continuous(expand = expansion(c(0, 0.04)),
                     sec.axis = sec_axis(~ . * WMAX / RMAX,
                                         name = expression(mean~-log[10]*italic(P)))) +
  coord_cartesian(xlim = c(0, xmax_all)) +
  labs(y = "cM / Mb") +
  theme_classic(base_size = 8) +
  theme(axis.title.x = element_blank(), legend.position = "none",
        axis.title.y.right = element_text(colour = "#D62728"),
        axis.text.y.right  = element_text(colour = "#D62728"),
        plot.margin = margin(6, 12, 2, 16))

png(OUT, width = W_IN, height = H_IN, units = "in", res = DPI)
print(p); dev.off()
cat("wrote", OUT, "\n")
