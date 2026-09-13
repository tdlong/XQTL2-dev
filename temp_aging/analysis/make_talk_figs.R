# make_talk_figs.R — talk versions of Figure 1 panels A and C.
#
#   Rscript temp_aging/analysis/make_talk_figs.R
#
# Writes temp_aging/analysis/Talk_Fig1.png (Wald, panel A) and Talk_Fig2.png
# (h2 partition, panel C). Same data and colours as make_figure1.R; the layout
# is five stacked rows, one per arm, so the plot is tall rather than wide.
#
# Talk_Fig3a.png and Talk_Fig3b.png are the allelic-effects slide, as two
# pieces to position separately: (a) a small genome-wide scan, split y axis, with
# the 2L and 3L peaks pinned; (b) the founder frequency changes under each. One
# treatment throughout (TRT3).
#
# Every row shares one x scale -- a cM is the same width in every row, so peak
# widths compare across arms, and arms shorter than X end short of the right
# edge. Every row shares one y scale too; only 3L, the one arm that leaves it,
# gets a broken axis.
#
# X_UNIT=Mb gives the physical axis instead of the default cM.

suppressMessages({library(tidyverse); library(patchwork)})

FOURSCAN <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
VARCOMP  <- "process/AGE_SY_splithalf/H2_varcomp_by_window_no89.txt.gz"
X_UNIT   <- Sys.getenv("X_UNIT", "cM")
W_IN <- 8; H_IN <- 7.5; DPI <- 300; BASE <- 14

CHRS   <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
CHRLAB <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")
HET <- read.table("pipeline/helpfiles/het_bounds.txt", header = TRUE,
                  comment.char = "#") %>% as_tibble()

TRT_LEV <- c("SY10 female", "SY20 female", "SY10 male", "SY20 male")
TRT_COL <- c("SY10 female" = "#F49AC2", "SY20 female" = "#D62728",
             "SY10 male"   = "#8EC7E8", "SY20 male"   = "#1F4E9C")
CMP_LEV <- c("main", "sex", "diet")
CMP_COL <- c(main = "grey25", sex = "#E69F00", diet = "#009E73")

# ── data ─────────────────────────────────────────────────────────────────────
scans <- read.table(FOURSCAN, header = TRUE, sep = "\t") %>% as_tibble() %>%
  filter(chr %in% CHRS) %>%
  transmute(chr, pos, cM, wald = Wald_log10p,
            trt = factor(paste0(sugar, " ", ifelse(sex == "F", "female", "male")),
                         levels = TRT_LEV))
posx <- function(pos, cM) if (X_UNIT == "cM") cM else pos / 1e6
scans <- scans %>% mutate(x = posx(pos, cM))

# the varcomp file has no cM: interpolate it per arm from the scan's pos -> cM
cmap <- scans %>% distinct(chr, pos, cM) %>% arrange(chr, pos)
vc <- read.table(VARCOMP, header = TRUE, sep = "\t") %>% as_tibble() %>%
  filter(chr %in% CHRS) %>% group_split(chr) %>% map_dfr(function(v) {
    m <- cmap %>% filter(chr == v$chr[1])
    v %>% mutate(cM = approx(m$pos, m$cM, xout = pos, rule = 2)$y)
  })
# panel C as in make_figure1.R: 100 kb bin means, no further smoothing
BIN <- 1e5
cmp <- vc %>% mutate(bin = pos %/% BIN) %>% group_by(chr, bin) %>%
  summarise(main = mean(mainH2), sex = mean(sexH2), diet = mean(dietH2),
            pos = mean(pos), cM = mean(cM), .groups = "drop") %>%
  pivot_longer(all_of(CMP_LEV), names_to = "term", values_to = "y") %>%
  mutate(term = factor(term, levels = CMP_LEV), x = posx(pos, cM))

# ── the shared x scale ───────────────────────────────────────────────────────
arms <- scans %>% group_by(chr) %>%
  summarise(lo = min(x), hi = max(x), .groups = "drop")
if (X_UNIT == "Mb") arms$lo <- 0
SPAN <- max(arms$hi - arms$lo)
het <- HET %>% filter(chr %in% CHRS) %>% pmap_dfr(function(chr, eu_start, eu_end) {
  m <- cmap %>% filter(chr == !!chr)
  at <- if (X_UNIT == "cM") approx(m$pos, m$cM, xout = c(eu_start, eu_end) * 1e6,
                                   rule = 2)$y else c(eu_start, eu_end)
  a <- arms %>% filter(chr == !!chr)
  tibble(chr, xmin = c(a$lo, at[2]), xmax = c(at[1], a$hi))
}) %>% filter(xmax > xmin)

# ── one row per arm; 3L split where it leaves the common scale ───────────────
# The break machinery is make_figure1.R's: sub-panel ranges exclusive and
# exhaustive (expand = 0), and data outside a segment set to NA, because the
# break marks need clip = "off" and that stops ggplot clipping the data.
seg <- function(d, yv, r) {
  d[[yv]][!is.na(d[[yv]]) & (d[[yv]] < r[1] | d[[yv]] > r[2])] <- NA_real_
  d
}
# Masking alone leaves a gap at the break: the last point below it and the first
# above it sit in different sub-panels, so neither draws the step between them.
# Add a point exactly at the break wherever a line crosses it; seg() keeps it in
# both sub-panels, so each line runs to the edge of both.
cross <- function(d, yv, cv, brk) {
  add <- d %>% arrange(.data[[cv]], x) %>% group_by(chr, .data[[cv]]) %>%
    mutate(y0 = .data[[yv]], y1 = lead(y0), x1 = lead(x)) %>% ungroup() %>%
    filter(!is.na(y1), (y0 - brk) * (y1 - brk) < 0) %>%
    transmute(chr, !!cv := .data[[cv]], x = x + (brk - y0) * (x1 - x) / (y1 - y0),
              !!yv := brk)
  bind_rows(d, add)
}

# one stroke where a sub-panel's y axis ends; `rel` is its height in the layout,
# so the stroke is the same size on the page in both sub-panels
slash <- function(yr, at, rel) {
  h <- diff(yr) * 0.05 / rel; w <- SPAN * 0.006
  annotate("segment", x = -w, xend = w, y = at - h, yend = at + h,
           colour = "black", linewidth = 0.5)
}

row_plot <- function(d, a, yv, cv, cols, yr, ybr, ylab, lw, zero,
                     top = FALSE, bottom_row = FALSE, rel = 1, brk_at = NULL) {
  lo <- arms$lo[arms$chr == a]
  dd <- seg(d %>% filter(chr == a), yv, yr) %>% mutate(x = x - lo)
  hb <- het %>% filter(chr == a) %>% mutate(xmin = xmin - lo, xmax = xmax - lo)
  # real coordinates on the ticks, even though each row is drawn from its own origin
  lab <- function(b) format(round(b + lo), trim = TRUE)
  show_x <- !top && (X_UNIT == "cM" || bottom_row)
  p <- ggplot() +
    geom_rect(data = hb, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey88", colour = NA) +
    { if (zero && yr[1] <= 0 && yr[2] >= 0)
        geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) } +
    geom_line(data = dd, aes(x, .data[[yv]], colour = .data[[cv]],
                             group = .data[[cv]]), linewidth = lw, na.rm = TRUE) +
    scale_colour_manual(values = cols, drop = FALSE, name = NULL) +
    scale_x_continuous(expand = expansion(0), labels = lab,
                       breaks = if (X_UNIT == "Mb") seq(0, 30, 5) else
                                  seq(-lo %% 10, SPAN, 10)) +
    scale_y_continuous(expand = expansion(0), breaks = ybr) +
    coord_cartesian(xlim = c(0, SPAN), ylim = yr, clip = "off") +
    labs(x = if (bottom_row) paste0("Position on arm (", X_UNIT, ")") else NULL,
         y = ylab, tag = if (top) NULL else CHRLAB[[a]]) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.5))) +
    theme_classic(base_size = BASE) +
    theme(plot.tag.position = "right",
          plot.tag = element_text(face = "bold", size = BASE),
          axis.text = element_text(size = BASE - 3),
          plot.margin = margin(if (top) 4 else 1, 4, if (top) 1 else 3, 4))
  if (!show_x) p <- p + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())
  if (top) p <- p + theme(axis.line.x = element_blank())
  if (!is.null(brk_at)) p <- p + slash(yr, brk_at, rel)
  p
}

stack <- function(d, yv, cv, cols, lo_r, lo_br, hi_r, hi_br, ylab, lw, zero,
                  HI_REL = 0.55) {
  d <- cross(d, yv, cv, lo_r[2])
  plots <- list(); hts <- c()
  for (a in CHRS) {
    last <- a == tail(CHRS, 1)
    if (a == "chr3L") {
      plots <- c(plots, list(row_plot(d, a, yv, cv, cols, hi_r, hi_br, ylab, lw, zero,
                                      top = TRUE, rel = HI_REL, brk_at = hi_r[1])))
      hts <- c(hts, HI_REL)
      plots <- c(plots, list(row_plot(d, a, yv, cv, cols, lo_r, lo_br, ylab, lw, zero,
                                      rel = 1, brk_at = lo_r[2])))
    } else {
      plots <- c(plots, list(row_plot(d, a, yv, cv, cols, lo_r, lo_br, ylab, lw, zero,
                                      bottom_row = last)))
    }
    hts <- c(hts, 1)
  }
  # the rows share one legend and one y title
  wrap_plots(plots, ncol = 1, heights = hts) +
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = "top", legend.text = element_text(size = BASE - 1),
          legend.key.width = unit(22, "pt"))
}

# what leaves the common scale must be 3L alone, or the layout is wrong
check <- function(d, yv, brk) {
  over <- d %>% filter(.data[[yv]] > brk) %>% distinct(chr) %>% pull(chr)
  if (!identical(sort(as.character(over)), "chr3L"))
    stop("arms above the break at ", brk, ": ", paste(over, collapse = ", "), call. = FALSE)
}

save <- function(p, f) {
  png(f, width = W_IN, height = H_IN, units = "in", res = DPI)
  print(p); invisible(dev.off()); cat("Wrote:", f, "\n")
}

check(scans, "wald", 45)
save(stack(scans, "wald", "trt", TRT_COL, c(0, 45), c(0, 20, 40),
           c(45, 170), c(100, 150), expression(-log[10] * italic(P)), 0.45, FALSE),
     "temp_aging/analysis/Talk_Fig1.png")

check(cmp, "y", 1.25)
save(stack(cmp, "y", "term", CMP_COL, c(-0.25, 1.25), c(0, 0.5, 1),
           c(1.25, 3.2), c(2, 3), expression(italic(h)^2 ~ "(%)"), 0.6, TRUE),
     "temp_aging/analysis/Talk_Fig2.png")

# ── Talk_Fig3: allelic effects at the two peaks ──────────────────────────────
MEANS <- "process/AGE_SY/AGE_SY_zoom_means.txt.gz"
TRT3  <- c(sugar = "SY20", sex = "M")      # both peaks are strong here: 27.2, 163.8
PEAKS <- tibble(n = 1:2, locus = c("chr2L:11.95", "chr3L:9.31"))
WIN3  <- 0.5e6                             # half-width drawn, as in Figure 2
RARE  <- 0.025                             # founders under this in controls are faded
FOUNDER_COL <- c(AB8 = "#E69F00", B1 = "#56B4E9", B2 = "#009E73", B3 = "#F0E442",
                 B4 = "#D55E00", B5 = "#CC79A7", B6 = "#000000", B7 = "#990099")
PIN <- "grey30"
trt_name <- paste(TRT3[["sugar"]], if (TRT3[["sex"]] == "M") "male" else "female")

means <- read.table(MEANS, header = TRUE, sep = "\t") %>% as_tibble() %>%
  filter(!REP %in% c(8, 9), sugar == TRT3[["sugar"]], sex == TRT3[["sex"]],
         locus %in% PEAKS$locus)
dfreq <- means %>%
  pivot_wider(names_from = TRT, values_from = freq, names_prefix = "f") %>%
  group_by(locus, chr, peak_pos, pos, founder) %>%
  summarise(D = mean(fZ - fC), fC = mean(fC), .groups = "drop") %>%
  filter(abs(pos - peak_pos) <= WIN3) %>%
  group_by(locus, founder) %>% mutate(rare = mean(fC) < RARE) %>% ungroup() %>%
  arrange(!rare)                           # faded lines drawn behind

# the scan for that treatment, arms end to end in Mb
g3 <- read.table(FOURSCAN, header = TRUE, sep = "\t") %>% as_tibble() %>%
  filter(chr %in% CHRS, sugar == TRT3[["sugar"]], sex == TRT3[["sex"]]) %>%
  mutate(chr = factor(chr, levels = CHRS))
off <- g3 %>% group_by(chr) %>% summarise(len = max(pos), .groups = "drop") %>%
  mutate(off = lag(cumsum(len + 1e6), default = 0))
g3 <- g3 %>% left_join(off, by = "chr") %>% mutate(gx = (pos + off) / 1e6)
pins <- PEAKS %>% left_join(dfreq %>% distinct(locus, chr, peak_pos), by = "locus") %>%
  mutate(chr = factor(chr, levels = CHRS)) %>% left_join(off, by = "chr") %>%
  mutate(gx = (peak_pos + off) / 1e6,
         top = map2_dbl(chr, peak_pos, ~ max(g3$Wald_log10p[g3$chr == .x &
                                                abs(g3$pos - .y) < 5e4])))
# Split y axis, as in Talk_Fig1: without it the 3L peak flattens 2L into the
# floor. Only 3L crosses the break in this treatment.
BRK3 <- 45; LO3 <- c(0, BRK3); HI3 <- c(BRK3, 215); REL3 <- c(hi = 1, lo = 2)
check(g3 %>% rename(wald = Wald_log10p), "wald", BRK3)
g3 <- g3 %>% mutate(x = gx, shade = factor(as.integer(chr) %% 2 == 0),
                    wald = Wald_log10p) %>% cross("wald", "shade", BRK3)
# pins sit a fixed distance above the peak ON THE PAGE, so the offset in data
# units depends on which segment the peak falls in
pins <- pins %>% mutate(seg = ifelse(top > BRK3, "hi", "lo"),
                        unit = ifelse(seg == "hi", diff(HI3) / REL3[["hi"]],
                                                   diff(LO3) / REL3[["lo"]]),
                        tri = top + 0.13 * unit, lab = top + 0.33 * unit)
SPAN3 <- max(g3$x)

manh_seg <- function(r, which) {
  pp <- pins %>% filter(seg == which)
  p <- ggplot(seg(g3, "wald", r), aes(x, wald, colour = shade, group = chr)) +
    geom_line(linewidth = 0.5, na.rm = TRUE) +
    scale_colour_manual(values = c(`TRUE` = "grey20", `FALSE` = "grey60"), guide = "none") +
    geom_point(data = pp, aes(gx, tri), inherit.aes = FALSE,
               shape = 25, size = 3.2, fill = PIN, colour = PIN) +
    geom_label(data = pp, aes(gx, lab, label = n), inherit.aes = FALSE,
               fill = PIN, colour = "white", fontface = "bold", size = 5,
               linewidth = 0, label.r = unit(2, "pt")) +
    scale_x_continuous(breaks = (off$off + off$len / 2) / 1e6, labels = CHRLAB[CHRS],
                       expand = expansion(0)) +
    scale_y_continuous(breaks = if (which == "hi") c(100, 150) else c(0, 20, 40),
                       expand = expansion(0)) +
    coord_cartesian(xlim = c(0, SPAN3), ylim = r, clip = "off") +
    labs(x = NULL, y = expression(-log[10] * italic(P))) +
    theme_classic(base_size = BASE) +
    theme(axis.ticks.x = element_blank(), axis.text = element_text(size = BASE - 2))
  h <- diff(r) * 0.05 / REL3[[which]]; w <- SPAN3 * 0.006
  at <- if (which == "hi") r[1] else r[2]
  p <- p + annotate("segment", x = -w, xend = w, y = at - h, yend = at + h,
                    colour = "black", linewidth = 0.5)
  # the treatment label goes in the upper segment, over X, which is empty there
  if (which == "hi") p + annotate("label", x = SPAN3 * 0.015, y = r[2], label = trt_name,
                                  hjust = 0, vjust = 1, fill = TRT_COL[[trt_name]],
                                  colour = "white", size = 5, linewidth = 0) +
    theme(axis.text.x = element_blank(), axis.line.x = element_blank(),
          plot.margin = margin(4, 4, 1, 4)) else
    p + theme(plot.margin = margin(1, 4, 3, 4))
}
manh <- manh_seg(HI3, "hi") / manh_seg(LO3, "lo") +
  plot_layout(heights = REL3[c("hi", "lo")], axis_titles = "collect")

freq3 <- function(k) {
  lc <- PEAKS$locus[k]; d <- dfreq %>% filter(locus == lc)
  pk <- d$peak_pos[1] / 1e6; yr <- range(d$D)
  ggplot(d, aes(pos / 1e6, D, colour = founder, alpha = ifelse(rare, 0.25, 1))) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = pk, linetype = "dashed", linewidth = 0.4) +
    geom_line(linewidth = 1) +
    annotate("label", x = pk, y = yr[2] + diff(yr) * 0.14, label = PEAKS$n[k],
             fill = PIN, colour = "white", fontface = "bold", size = 5,
             linewidth = 0, label.r = unit(2, "pt")) +
    scale_colour_manual(values = FOUNDER_COL, name = NULL) + scale_alpha_identity() +
    scale_x_continuous(expand = expansion(0)) +
    scale_y_continuous(limits = c(yr[1] - diff(yr) * 0.04, yr[2] + diff(yr) * 0.2)) +
    labs(x = paste0(CHRLAB[[d$chr[1]]], " position (Mb)"),
         y = "Frequency change\n(aged − control)") +
    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 2))) +
    theme_classic(base_size = BASE) +
    theme(axis.text = element_text(size = BASE - 2))
}

save3 <- function(p, f, w, h) {
  png(f, width = w, height = h, units = "in", res = DPI)
  print(p); invisible(dev.off()); cat("Wrote:", f, "\n")
}
save3(manh, "temp_aging/analysis/Talk_Fig3a.png", 7, 3)
save3(freq3(1) + freq3(2) +
        plot_layout(axis_titles = "collect", guides = "collect") &
        theme(legend.position = "bottom", legend.text = element_text(size = BASE - 1),
              legend.key.width = unit(22, "pt")),
      "temp_aging/analysis/Talk_Fig3b.png", 10, 4.5)
