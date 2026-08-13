# make_figure1.R — Figure 1: Wald scan, Cutler h2, and the h2 partition.
#
# Three stacked Manhattan panels sharing one concatenated x-axis.
#   A  Wald -log10(p), four treatments
#   B  Cutler h2, four treatments
#   C  h2 partitioned into main, sex and diet
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript temp_aging/make_figure1.R
#
# SCANS: panels A and B are drawn from the split-half runs, averaged over the
# odd and even halves, because that is what is held locally. Those are 6
# replicates per half, so the Wald peaks are lower than the 12-replicate scans
# behind AGE_SY_v6_size75k_4scan.png. To use the 12-replicate scans instead,
# point LIVE_SCAN_DIR at a folder holding <scan>/<scan>.scan.txt for the four
# treatments; everything else is unchanged.

suppressMessages({library(tidyverse); library(patchwork)})

LONG      <- "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
VARCOMP   <- "process/AGE_SY_splithalf/H2_varcomp_by_window.txt.gz"
# The 12-replicate scans, collapsed to one small file by make_4scan_df.R on the
# cluster. If it is present it is used for panels A and B; otherwise those fall
# back to the split-half scans averaged over halves, which is NOT the same thing
# (6 replicates per half, so the Wald peaks are lower).
FOURSCAN <- "process/AGE_SY/AGE_SY_4scan.txt.gz"
OUT       <- "temp_aging/Figure1_plot.png"
W_IN <- 7.5; H_IN <- 6; DPI <- 300
SMOOTH_BP_C <- 5e5             # rolling mean on panel C only
YLIM_C <- c(-0.1, 1)

CHRS   <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
CHRLAB <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")
# dm6 euchromatin boundaries (Huynh et al. 2023 PLoS Genet 19:e1010439, Table S2)
HET <- tribble(~chr, ~eu_start, ~eu_end,
  "chrX", 2.5, 21.2, "chr2L", 0.5, 22.9, "chr2R", 1.3, 25.1,
  "chr3L", 0.7, 24.0, "chr3R", 4.5, 32.0)

TRT_LEV <- c("SY10 female", "SY20 female", "SY10 male", "SY20 male")
TRT_COL <- c("SY10 female" = "#F49AC2", "SY20 female" = "#D62728",
             "SY10 male"   = "#8EC7E8", "SY20 male"   = "#1F4E9C")
CMP_LEV <- c("main", "sex", "diet")
CMP_COL <- c(main = "grey25", sex = "#E69F00", diet = "#009E73")

# ── data ─────────────────────────────────────────────────────────────────────
long <- read.table(LONG, header = TRUE, sep = "\t") %>% as_tibble()
vc   <- read.table(VARCOMP, header = TRUE, sep = "\t") %>% as_tibble()

scans <- if (file.exists(FOURSCAN)) {
  cat("panels A/B: 12-replicate scans from", FOURSCAN, "\n")
  read.table(FOURSCAN, header = TRUE, sep = "\t") %>% as_tibble() %>%
    transmute(chr, pos, sugar, sex, wald = Wald_log10p, h2 = Cutl_H2)
} else {
  cat("panels A/B: FALLING BACK to the split-half scans averaged over halves.\n",
      "  Run temp_aging/make_4scan_df.R on hpc3 and fetch", basename(FOURSCAN), "\n", sep = "")
  long %>% group_by(chr, pos, sugar, sex) %>%
    summarise(wald = mean(Wald_log10p), h2 = mean(Cutl_H2), .groups = "drop")
}
scans <- scans %>%
  mutate(trt = factor(paste0(sugar, " ", ifelse(sex == "F", "female", "male")),
                      levels = TRT_LEV))

# ── one concatenated x-axis shared by all three panels ───────────────────────
lens <- scans %>% filter(chr %in% CHRS) %>% group_by(chr) %>%
  summarise(len = max(pos), .groups = "drop") %>%
  mutate(chr = factor(chr, levels = CHRS)) %>% arrange(chr) %>%
  mutate(offset = lag(cumsum(len), default = 0))

addx <- function(d) d %>% filter(chr %in% CHRS) %>%
  mutate(chr = factor(chr, levels = CHRS)) %>%
  left_join(lens %>% select(chr, offset), by = "chr") %>%
  mutate(gx = (pos + offset) / 1e6)

scans <- addx(scans)
vcx   <- addx(vc)

het_bands <- HET %>% mutate(chr = factor(chr, levels = CHRS)) %>%
  left_join(lens, by = "chr") %>%
  { bind_rows(
      transmute(., xmin = offset/1e6, xmax = (offset + eu_start*1e6)/1e6),
      transmute(., xmin = (offset + eu_end*1e6)/1e6, xmax = (offset + len)/1e6)) } %>%
  filter(xmax > xmin)
chr_breaks <- (lens$offset + lens$len/2) / 1e6           # label at chromosome centre
chr_edges  <- (lens$offset + lens$len)[-nrow(lens)] / 1e6 # divider between chromosomes
xmax_all   <- sum(lens$len) / 1e6

# ── panel C components, lightly smoothed ─────────────────────────────────────
roll <- function(x, k) { n <- length(x); h <- (k-1L) %/% 2L
  vapply(seq_len(n), function(i) mean(x[max(1L,i-h):min(n,i+h)], na.rm=TRUE), numeric(1)) }
BIN <- 1e5
cmp <- vcx %>%
  mutate(bin = (pos %/% BIN) * BIN) %>%
  group_by(chr, offset, bin) %>%
  summarise(main = mean(mainH2), sex = mean(sexH2), diet = mean(dietH2), .groups = "drop") %>%
  pivot_longer(all_of(CMP_LEV), names_to = "term", values_to = "y") %>%
  arrange(chr, term, bin) %>% group_by(chr, term) %>%
  mutate(y = roll(y, max(1L, as.integer(round(SMOOTH_BP_C / BIN))))) %>% ungroup() %>%
  mutate(term = factor(term, levels = CMP_LEV), gx = (bin + offset) / 1e6)

# ── panels ───────────────────────────────────────────────────────────────────
# Broken y axis: each panel is two ggplots stacked 3:1, so the chr3L locus does
# not set the scale for the whole genome.
#
# Not ggbreak. It renders correctly on its own but cannot be composed: through
# aplot::plot_list() the sub-break data is drawn again, squashed along the
# bottom of the upper segment (the traces on 3L/3R in panel A), and through
# ggplotify + cowplot the break is silently dropped altogether.
#
# Two things are needed and both matter:
#
#   expand = expansion(0)  -- the sub-panel ranges must be mutually exclusive
#     and exhaustive. With the default 5% expansion c(0, 20) actually renders as
#     -1 to 21 and c(20, 105) as ~15.7 to 109, so they overlap and anything in
#     the overlap is drawn in BOTH sub-panels.
#
#   mask()                 -- the break marks need clip = "off", which also
#     stops out-of-range data being clipped, so it spills across the whole
#     figure. Setting the wrong side of the break to NA keeps geom_line from
#     drawing it at all.

mask <- function(d, yv, keep_above, brk) {
  d[[yv]][if (keep_above) d[[yv]] < brk else d[[yv]] > brk] <- NA_real_
  d
}

deco <- function(xaxis = FALSE) {
  th <- theme_classic(base_size = 8) +
    theme(axis.title.x = element_blank(),
          plot.tag = element_text(size = 10, face = "bold"),
          plot.tag.position = c(0.004, 0.88),
          legend.position = "none")
  if (!xaxis) th <- th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  list(
    geom_rect(data = het_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey88", colour = NA, inherit.aes = FALSE),
    geom_vline(xintercept = chr_edges, linetype = "dotted", colour = "grey35", linewidth = 0.3),
    scale_x_continuous(expand = expansion(0),
                       breaks = chr_breaks, labels = CHRLAB[CHRS]),
    th)
}

# marks the break with two short diagonals on the y axis
# the // that marks the discontinuity, straddling the y axis
# ONE stroke, marking where this sub-panel's y axis terminates. The lower panel
# gets one at its top and the upper panel one at its bottom; with the gap
# between the panels those two strokes are the break notation. (Two strokes per
# panel gives four at the break, which is wrong.)
# `rel` is the sub-panel's relative height in the layout (1 for the short upper
# segment, 3 for the tall lower one). Dividing by it makes the stroke the same
# size ON THE PAGE in both, hence the same angle -- a fixed fraction of each
# panel's data range would make the lower stroke three times taller.
slash <- function(yr, at, rel) {
  h <- diff(yr) * 0.030 / rel
  w <- xmax_all * 0.004
  annotate("segment", x = -w, xend = w, y = at - h, yend = at + h,
           colour = "black", linewidth = 0.4)
}

split_panel <- function(d, yv, cv, cols, brk, lo, hi, lo_brk, hi_brk,
                        ylab, tag, lw = 0.3, zero = FALSE, xaxis = FALSE) {
  # group by chromosome as well as series, or geom_line joins the last window of
  # one arm to the first of the next straight through the centromere
  gl <- function(dd) geom_line(data = dd,
                               aes(gx, .data[[yv]], colour = .data[[cv]],
                                   group = interaction(.data[[cv]], chr)),
                               linewidth = lw)
  # the zero line belongs only to the segment whose range contains zero; with
  # clip = "off" the other segment would draw it far outside its own panel
  zl <- function(r) if (zero && r[1] <= 0 && r[2] >= 0)
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) else NULL
  top <- ggplot() + deco() + zl(hi) + gl(mask(d, yv, TRUE, brk)) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_y_continuous(breaks = hi_brk, expand = expansion(0)) +
    coord_cartesian(xlim = c(0, xmax_all), ylim = hi, clip = "off") + slash(hi, hi[1], 1) +
    labs(y = NULL, tag = tag) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), plot.margin = margin(7, 4, 2, 16))
  bot <- ggplot() + deco(xaxis) + zl(lo) + gl(mask(d, yv, FALSE, brk)) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_y_continuous(breaks = lo_brk, expand = expansion(0)) +
    coord_cartesian(xlim = c(0, xmax_all), ylim = lo, clip = "off") + slash(lo, lo[2], 3) +
    labs(y = ylab) +
    theme(plot.margin = margin(2, 4, 2, 16))
  list(top = top, bot = bot)
}

pA <- split_panel(scans, "wald", "trt", TRT_COL, 20, c(0, 20), c(20, 215),
                  c(0, 5, 10, 15), c(50, 100, 150, 200),
                  expression(-log[10] * italic(P)), "A", lw = 0.28)
pB <- split_panel(scans, "h2", "trt", TRT_COL, 2.5, c(0, 2.5), c(2.5, 5),
                  c(0, 0.5, 1, 1.5, 2), c(3, 4, 5),
                  expression(italic(h)^2), "B", lw = 0.28)
pC <- split_panel(cmp, "y", "term", CMP_COL, 0.75, c(-0.25, 0.75), c(0.75, 2.9),
                  c(-0.2, 0, 0.25, 0.5), c(1, 2),
                  expression(italic(h)^2), "C", lw = 0.35, zero = TRUE, xaxis = TRUE)

key <- function(lev, col) {
  ggplot(tibble(x = 1, y = 1, g = factor(lev, levels = lev)), aes(x, y, colour = g)) +
    geom_line() + scale_colour_manual(values = col, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(nrow = 1, keywidth = unit(9, "pt"),
                                 override.aes = list(linewidth = 0.8))) +
    theme_void(base_size = 8) + theme(legend.position = "bottom")
}

# one patchwork over all six sub-panels: patchwork aligns the panel areas, so
# the three groups line up even though their y labels differ in width
fig <- pA$top / pA$bot / pB$top / pB$bot / pC$top / pC$bot +
  plot_layout(heights = c(1, 3, 1, 3, 1, 3))

g <- cowplot::plot_grid(
  fig,
  cowplot::plot_grid(cowplot::get_legend(key(TRT_LEV, TRT_COL)),
                     cowplot::get_legend(key(CMP_LEV, CMP_COL)),
                     nrow = 1, rel_widths = c(2.1, 1)),
  ncol = 1, rel_heights = c(1, 0.075))

png(OUT, width = W_IN, height = H_IN, units = "in", res = DPI)
print(g)
invisible(dev.off())
cat("Wrote:", OUT, "\n")
