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
#
# DATA: the 10-replicate AGE_SY dataset -- replicates 8 and 9 dropped, those
# being the May 2023 cage (helpfiles/AGE_2024/population_assignment.txt). The
# 12-replicate files still exist on disk; nothing here reads them.

suppressMessages({library(tidyverse); library(patchwork)})
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# Optional variant, e.g.  Rscript temp_aging/make_figure1.R no89
# "no89" is AGE_SY without replicates 8 and 9 -- the May 2023 cage -- leaving a
# single-cage experiment on 10 replicates that still halves evenly 5/5. Inputs
# gain the suffix, the figure is written as Figure1b.

LONG      <- "process/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz"
VARCOMP   <- "process/AGE_SY_splithalf/H2_varcomp_by_window_no89.txt.gz"
# The 12-replicate scans, collapsed to one small file by make_4scan_df.R on the
# cluster. If it is present it is used for panels A and B; otherwise those fall
# back to the split-half scans averaged over halves, which is NOT the same thing
# (6 replicates per half, so the Wald peaks are lower).
FOURSCAN <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
# X_UNIT "Mb" gives the physical axis; X_UNIT=cM gives the genetic one, written
# as Figure1_cM_plot.png. In cM, 2L and 2R are contiguous (2L ends at 54.0, 2R
# starts at 54.6) and so are 3L and 3R, so that axis concatenates by LINKAGE
# GROUP -- X, 2, 3 -- with the arm boundaries still marked inside each. Peak
# width in cM is the thing that should be roughly constant; in Mb it tracks
# recombination.
X_UNIT    <- Sys.getenv("X_UNIT", "Mb")
OUT       <- sprintf("temp_aging/Figure1%s%s%s_plot.png",
                     if (toupper(Sys.getenv("PANELS","ABC"))=="AB") "_AB" else "",
                     if (Sys.getenv("H2","cutler")=="falconer") "_falc" else "",
                     if (X_UNIT == "cM") "_cM" else "")
W_IN <- 7.5; H_IN <- 6; DPI <- 300
SMOOTH_BP_C <- 5e5             # rolling mean on panel C only
YLIM_C <- c(-0.1, 1)

CHRS   <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
CHRLAB <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")
# Euchromatin boundaries: read from the pipeline rather than hardcoded. These
# define the grey bands in the figures; the analysis scripts also use them to
# restrict to euchromatin.
HET <- read.table("pipeline/helpfiles/het_bounds.txt", header = TRUE,
                  comment.char = "#") %>% as_tibble()

TRT_LEV <- c("SY10 female", "SY20 female", "SY10 male", "SY20 male")
TRT_COL <- c("SY10 female" = "#F49AC2", "SY20 female" = "#D62728",
             "SY10 male"   = "#8EC7E8", "SY20 male"   = "#1F4E9C")
CMP_LEV <- c("main", "sex", "diet")
CMP_COL <- c(main = "grey25", sex = "#E69F00", diet = "#009E73")

# ── data ─────────────────────────────────────────────────────────────────────
for (f in c(LONG, VARCOMP)) if (!file.exists(f)) stop("missing ", f,
  "
Run scripts_oneoffs/AGE_SY/nov_only/{run_scans.sh,gather.R} on HPC3, scp back,",
  "
then varcomp_H2.R with the two paths as arguments.", call. = FALSE)

# PANELS=AB draws the scan and the heritability only. Panel C is the variance
# partition, which needs the eight split-half scans and the varcomp file derived
# from them -- neither of which says anything about whether an h2 estimator is
# right, and both of which go stale the moment the main scans are rerun. While a
# fix is being validated, AB is the whole test.
PANELS <- toupper(Sys.getenv("PANELS", "ABC"))
want_C <- grepl("C", PANELS, fixed = TRUE)

long <- if (want_C || !file.exists(FOURSCAN))
  read.table(LONG, header = TRUE, sep = "\t") %>% as_tibble() else NULL
vc   <- if (want_C)
  read.table(VARCOMP, header = TRUE, sep = "\t") %>% as_tibble() else NULL

scans <- if (file.exists(FOURSCAN)) {
  cat("panels A/B:", basename(FOURSCAN), "\n")
  read.table(FOURSCAN, header = TRUE, sep = "\t") %>% as_tibble() %>%
    transmute(chr, pos, cM, sugar, sex, wald = Wald_log10p, h2 = Cutl_H2)
} else {
  cat("panels A/B: FALLING BACK to the split-half scans averaged over halves.\n",
      "  Run temp_aging/make_4scan_df.R on hpc3 and fetch", basename(FOURSCAN), "\n", sep = "")
  # the split-half file has no cM, so this fallback cannot draw the genetic axis
  long %>% group_by(chr, pos, sugar, sex) %>%
    summarise(wald = mean(Wald_log10p), h2 = mean(Cutl_H2),
              cM = NA_real_, .groups = "drop")
}
# Panel B source. Default is the pipeline's Cutler h2, uncorrected. H2=falconer
# swaps in the bias-corrected Falconer h2 from h2_from_scan.R (XQTL2 #40): the
# Cutler bias saturates against the penetrance clamp and so under-corrects where
# the variance is largest, which is the male X.
if (Sys.getenv("H2", "cutler") == "falconer") {
  fal <- map_dfr(c("SY10_F","SY20_F","SY10_M","SY20_M"), function(sc) {
    f <- file.path("process/AGE_SY/Scans", paste0("AGE_", sc, "_no89"),
                   paste0("AGE_", sc, "_no89.h2_falconer.txt"))
    if (!file.exists(f)) stop("missing ", f, "\nRun reproduce.sh on HPC3, then make_figures.sh.")
    read.table(f, header = TRUE, colClasses = c(sex = "character")) %>% as_tibble() %>%
      transmute(chr = CHROM, pos, sugar = str_extract(sc, "SY[12]0"),
                sex = str_sub(sc, -1), h2_new = h2_vc, bias_new = h2_bias)
  })
  # The Falconer bias carries a 1/C_f term, so where a founder sits at the lsei
  # floor it diverges -- p99 is 983 on chrX males, against an h2 scale where the
  # strongest locus in the genome is 4.3. Subtracting that gives h2_corr of -40,
  # which does not just look wrong, it destroys the panel scales and bleeds
  # outside the plot. Those windows carry no information about h2 either way, so
  # they are masked to NA and the line breaks there (XQTL2 #40).
  # h2_vc, not h2_corr: it fits tau2 per window by ML, so non-negativity is a
  # boundary solution rather than a clamp. h2_corr is unbiased but unbounded
  # below and came out negative at 85% of chrX windows (XQTL2 #40).
  scans <- scans %>% inner_join(fal, by = c("chr", "pos", "sugar", "sex")) %>%
    mutate(h2 = h2_new)
  cat(sprintf("panel B: variance-component h2 (h2_vc), %d windows; %.1f%% at zero\n",
              nrow(scans), 100*mean(scans$h2 <= 0, na.rm = TRUE)))
  scans <- scans %>% select(-h2_new, -bias_new)
} else {
  cat("panel B: Cutler h2, uncorrected\n")
}

scans <- scans %>%
  mutate(trt = factor(paste0(sugar, " ", ifelse(sex == "F", "female", "male")),
                      levels = TRT_LEV))

# ── one concatenated x-axis shared by all three panels ───────────────────────
if (X_UNIT == "cM") {
  if (all(is.na(scans$cM))) stop("no cM in the scan file", call. = FALSE)
  GROUP <- c(chrX = "X", chr2L = "2", chr2R = "2", chr3L = "3", chr3R = "3")
  GAP   <- 4                                # cM of white space between groups

  # within a linkage group cM is already continuous across the two arms, so the
  # offset is per group, not per arm
  grp <- scans %>% filter(chr %in% CHRS) %>% mutate(g = GROUP[as.character(chr)]) %>%
    group_by(g) %>% summarise(gmin = min(cM, na.rm = TRUE),
                              gmax = max(cM, na.rm = TRUE), .groups = "drop") %>%
    mutate(g = factor(g, levels = c("X", "2", "3"))) %>% arrange(g) %>%
    mutate(span = gmax - gmin,
           offset = lag(cumsum(span + GAP), default = 0) - gmin)

  addx <- function(d) d %>% filter(chr %in% CHRS) %>%
    mutate(chr = factor(chr, levels = CHRS), g = factor(GROUP[as.character(chr)],
                                                        levels = c("X","2","3"))) %>%
    left_join(grp %>% select(g, offset), by = "g") %>%
    mutate(gx = cM + offset)

  # vc carries no cM, so interpolate it from the scan's own pos -> cM per arm
  if (want_C) vc <- vc %>% group_split(chr) %>% map_dfr(function(v) {
    x <- scans %>% filter(chr == v$chr[1]) %>% arrange(pos)
    if (!nrow(x)) return(v %>% mutate(cM = NA_real_))
    v %>% mutate(cM = approx(x$pos, x$cM, xout = pos, rule = 2)$y)
  })
  scans <- addx(scans); vcx <- if (want_C) addx(vc) else NULL

  # arm extents in cM, taken from the scan itself so they cannot disagree with it
  arms <- scans %>% group_by(chr) %>%
    summarise(lo = min(gx), hi = max(gx), .groups = "drop") %>%
    mutate(chr = factor(chr, levels = CHRS)) %>% arrange(chr)

  # heterochromatin: interpolate each arm's euchromatin boundaries onto cM
  het_bands <- HET %>% pmap_dfr(function(chr, eu_start, eu_end) {
    x <- scans %>% filter(chr == !!chr) %>% arrange(pos)
    at <- approx(x$pos, x$gx, xout = c(eu_start, eu_end) * 1e6, rule = 2)$y
    tibble(chr, xmin = c(min(x$gx), at[2]), xmax = c(at[1], max(x$gx)))
  }) %>% filter(xmax > xmin)

  chr_breaks <- (arms$lo + arms$hi) / 2                # label each arm
  chr_edges  <- sort(c(arms$hi[-nrow(arms)]))          # arm and group boundaries
  xmax_all   <- max(scans$gx)
  XLAB       <- "cM"
} else {
  lens <- scans %>% filter(chr %in% CHRS) %>% group_by(chr) %>%
    summarise(len = max(pos), .groups = "drop") %>%
    mutate(chr = factor(chr, levels = CHRS)) %>% arrange(chr) %>%
    mutate(offset = lag(cumsum(len), default = 0))

  addx <- function(d) d %>% filter(chr %in% CHRS) %>%
    mutate(chr = factor(chr, levels = CHRS)) %>%
    left_join(lens %>% select(chr, offset), by = "chr") %>%
    mutate(gx = (pos + offset) / 1e6)

  scans <- addx(scans); vcx <- if (want_C) addx(vc) else NULL

  het_bands <- HET %>% mutate(chr = factor(chr, levels = CHRS)) %>%
    left_join(lens, by = "chr") %>%
    { bind_rows(
        transmute(., xmin = offset/1e6, xmax = (offset + eu_start*1e6)/1e6),
        transmute(., xmin = (offset + eu_end*1e6)/1e6, xmax = (offset + len)/1e6)) } %>%
    filter(xmax > xmin)
  chr_breaks <- (lens$offset + lens$len/2) / 1e6
  chr_edges  <- (lens$offset + lens$len)[-nrow(lens)] / 1e6
  xmax_all   <- sum(lens$len) / 1e6
  XLAB       <- "Mb"
}

# ── panel C components, lightly smoothed ─────────────────────────────────────
roll <- function(x, k) { n <- length(x); h <- (k-1L) %/% 2L
  vapply(seq_len(n), function(i) mean(x[max(1L,i-h):min(n,i+h)], na.rm=TRUE), numeric(1)) }
BIN <- 1e5
cmp <- if (!want_C) NULL else vcx %>%
  mutate(bin = (pos %/% BIN) * BIN) %>%
  group_by(chr, offset, bin) %>%
  summarise(main = mean(mainH2), sex = mean(sexH2), diet = mean(dietH2),
            bin_gx = mean(gx), .groups = "drop") %>%
  pivot_longer(all_of(CMP_LEV), names_to = "term", values_to = "y") %>%
  arrange(chr, term, bin) %>% group_by(chr, term) %>%
  mutate(y = roll(y, max(1L, as.integer(round(SMOOTH_BP_C / BIN))))) %>% ungroup() %>%
  mutate(term = factor(term, levels = CMP_LEV),
         gx = if (X_UNIT == "cM") bin_gx else (bin + offset) / 1e6)

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

# clip = "off" is needed for the break marks, and it also stops ggplot clipping
# out-of-range DATA -- a value below a segment's floor is drawn outside the
# panel, across the figure and into the legend. mask() only splits at the break,
# so it does not catch that. NA out anything outside the segment's own range.
clip_to <- function(d, yv, r) {
  d[[yv]][!is.na(d[[yv]]) & (d[[yv]] < r[1] | d[[yv]] > r[2])] <- NA_real_
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
  top <- ggplot() + deco() + zl(hi) + gl(clip_to(mask(d, yv, TRUE, brk), yv, hi)) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_y_continuous(breaks = hi_brk, expand = expansion(0)) +
    coord_cartesian(xlim = c(0, xmax_all), ylim = hi, clip = "off") + slash(hi, hi[1], 1) +
    labs(y = NULL, tag = tag) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), plot.margin = margin(7, 4, 2, 16))
  bot <- ggplot() + deco(xaxis) + zl(lo) + gl(clip_to(mask(d, yv, FALSE, brk), yv, lo)) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_y_continuous(breaks = lo_brk, expand = expansion(0)) +
    coord_cartesian(xlim = c(0, xmax_all), ylim = lo, clip = "off") + slash(lo, lo[2], 3) +
    labs(y = ylab) +
    theme(plot.margin = margin(2, 4, 2, 16))
  list(top = top, bot = bot)
}

pA <- split_panel(scans, "wald", "trt", TRT_COL, 45, c(0, 45), c(45, 215),
                  c(0, 10, 20, 30, 40), c(50, 100, 150, 200),
                  expression(-log[10] * italic(P)), "A", lw = 0.28)
# The corrected h2 is an unbiased estimate of zero at null windows, so it
# scatters negative and the panel has to show that. The Cutler version cannot go
# below zero, hence the floor of 0 in the default range.
# h2_vc is non-negative by construction, so the panel keeps its zero floor.
pB <- split_panel(scans, "h2", "trt", TRT_COL, 2.5, c(0, 2.5), c(2.5, 5),
                  c(0, 0.5, 1, 1.5, 2), c(3, 4, 5),
                  expression(italic(h)^2), "B", lw = 0.28, xaxis = !want_C)
# Break at 1.25, not 0.75: the pericentromeric main plateau runs to ~1.0, and a
# break below it pushed that whole stretch into the short upper sub-panel. Only
# 0.3-0.4% of smoothed values exceed 1.25 -- essentially just the chr3L spike --
# in both the 12- and 10-replicate versions.
pC <- if (!want_C) NULL else
  split_panel(cmp, "y", "term", CMP_COL, 1.25, c(-0.25, 1.25), c(1.25, 2.9),
                  c(-0.2, 0, 0.5, 1.0), c(1.5, 2, 2.5),
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
fig <- if (want_C)
  pA$top / pA$bot / pB$top / pB$bot / pC$top / pC$bot +
    plot_layout(heights = c(1, 3, 1, 3, 1, 3)) else
  pA$top / pA$bot / pB$top / pB$bot +
    plot_layout(heights = c(1, 3, 1, 3))

g <- cowplot::plot_grid(
  fig,
  cowplot::plot_grid(cowplot::get_legend(key(TRT_LEV, TRT_COL)),
                     if (want_C) cowplot::get_legend(key(CMP_LEV, CMP_COL)) else
                       ggplot() + theme_void(),
                     nrow = 1, rel_widths = c(2.1, 1)),
  ncol = 1, rel_heights = c(1, 0.075))

png(OUT, width = W_IN, height = H_IN, units = "in", res = DPI)
print(g)
invisible(dev.off())
cat("Wrote:", OUT, "\n")
