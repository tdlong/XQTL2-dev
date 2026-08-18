# make_figureX.R — the Aug 2024 pilot against AGE_SY, at matched replicate count.
#
# Two panels, built like Figure 1's A and B and sharing its colours:
#   A  Wald -log10(p)
#   B  Cutler h2
#
# Five series: the four AGE_SY treatments cut to 6 replicates, in the Figure 1
# colours, plus AGE_2024 (6 cages, females, lab food) in cyan.
#
# Writes two files:
#   FigureX_plot.png    all five
#   FigureXb_plot.png   the three females only -- three foods, one sex, which is
#                       the comparison the males only clutter
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript temp_aging/make_figureX.R
#
# Needs process/AGE_2024/AGE_2024_vs_AGE_SY.txt.gz, written by
# scripts_oneoffs/AGE_2024/gather_scans.R on HPC3 and fetched with scp.

suppressMessages({library(tidyverse); library(patchwork)})

GATHER <- "process/AGE_2024/AGE_2024_vs_AGE_SY.txt.gz"
OUT_A  <- "temp_aging/FigureX_plot.png"
OUT_B  <- "temp_aging/FigureXb_plot.png"
W_IN <- 7.5; H_IN <- 4.2; DPI <- 300

# Which 6-replicate cut of AGE_SY to draw. "1-6" is the contiguous block; "odd"
# and "even" are the split-half scans, and switching to them is a one-word edit.
SY_REPS <- "1-6"

# y-axis breaks: lower range, upper range, and the ticks for each. Set BRK to NA
# to draw an unbroken panel instead. The 6-replicate peaks are lower than the
# 12-replicate ones behind Figure 1, so these are not Figure 1's numbers -- the
# script prints the observed maxima so they can be checked.
BRK_A <- 25; LO_A <- c(0, 25);  HI_A <- c(25, 105)
BRK_B <- 2.5; LO_B <- c(0, 2.5); HI_B <- c(2.5, 5.5)
TICK_LO_A <- c(0, 5, 10, 15, 20); TICK_HI_A <- c(40, 60, 80, 100)
TICK_LO_B <- c(0, 0.5, 1, 1.5, 2); TICK_HI_B <- c(3, 4, 5)

CHRS   <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
CHRLAB <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")
# dm6 euchromatin boundaries (Huynh et al. 2023 PLoS Genet 19:e1010439, Table S2)
HET <- tribble(~chr, ~eu_start, ~eu_end,
  "chrX", 2.5, 21.2, "chr2L", 0.5, 22.9, "chr2R", 1.3, 25.1,
  "chr3L", 0.7, 24.0, "chr3R", 4.5, 32.0)

# Figure 1's four, plus the pilot. The palette's logic is warm = female, cool =
# male, light = SY10, dark = SY20, so the third female is a warm colour: cyan
# read as a male, sitting right next to SY10 male's light blue.
LEV <- c("SY10 female", "SY20 female", "lab female", "SY10 male", "SY20 male")
COL <- c("SY10 female" = "#F49AC2", "SY20 female" = "#D62728",
         "lab female"  = "#FF7F0E",
         "SY10 male"   = "#8EC7E8", "SY20 male"   = "#1F4E9C")
FEMALES <- c("SY10 female", "SY20 female", "lab female")

# ── data ─────────────────────────────────────────────────────────────────────
if (!file.exists(GATHER))
  stop("missing ", GATHER,
       "\nRun scripts_oneoffs/AGE_2024/gather_scans.R on HPC3 and scp it back.")

g <- read.table(GATHER, header = TRUE, sep = "\t") %>% as_tibble()

scans <- g %>%
  filter((diet != "lab" & reps == SY_REPS) | diet == "lab") %>%
  mutate(trt = factor(case_when(
           diet == "lab" ~ "lab female",
           TRUE ~ paste0(diet, " ", ifelse(sex == "F", "female", "male"))),
         levels = LEV)) %>%
  filter(!is.na(trt)) %>%
  transmute(chr, pos, trt, wald = Wald_log10p, h2 = Cutl_H2)

got <- scans %>% count(trt, name = "windows")
cat("series drawn (AGE_SY at reps '", SY_REPS, "'):\n", sep = "")
print(as.data.frame(got), row.names = FALSE)
if (nrow(got) < 5) cat("\nNOTE: expected 5 series -- some scans are missing.\n\n")

cat("\nobserved maxima -- check the y-axis breaks against these:\n")
scans %>% group_by(trt) %>%
  summarise(max_wald = round(max(wald, na.rm = TRUE), 1),
            max_h2 = round(max(h2, na.rm = TRUE), 2), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n")

# ── one concatenated x-axis, shared by both panels ───────────────────────────
lens <- scans %>% filter(chr %in% CHRS) %>% group_by(chr) %>%
  summarise(len = max(pos), .groups = "drop") %>%
  mutate(chr = factor(chr, levels = CHRS)) %>% arrange(chr) %>%
  mutate(offset = lag(cumsum(len), default = 0))

scans <- scans %>% filter(chr %in% CHRS) %>%
  mutate(chr = factor(chr, levels = CHRS)) %>%
  left_join(lens %>% select(chr, offset), by = "chr") %>%
  mutate(gx = (pos + offset) / 1e6)

het_bands <- HET %>% mutate(chr = factor(chr, levels = CHRS)) %>%
  left_join(lens, by = "chr") %>%
  { bind_rows(
      transmute(., xmin = offset/1e6, xmax = (offset + eu_start*1e6)/1e6),
      transmute(., xmin = (offset + eu_end*1e6)/1e6, xmax = (offset + len)/1e6)) } %>%
  filter(xmax > xmin)
chr_breaks <- (lens$offset + lens$len/2) / 1e6
chr_edges  <- (lens$offset + lens$len)[-nrow(lens)] / 1e6
xmax_all   <- sum(lens$len) / 1e6

# ── panels (same construction as Figure 1) ───────────────────────────────────
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
    scale_x_continuous(expand = expansion(0), breaks = chr_breaks, labels = CHRLAB[CHRS]),
    th)
}

# one stroke marking where a sub-panel's y axis terminates; `rel` is its relative
# height so the stroke is the same size on the page in both, hence same angle
slash <- function(yr, at, rel) {
  h <- diff(yr) * 0.030 / rel
  w <- xmax_all * 0.004
  annotate("segment", x = -w, xend = w, y = at - h, yend = at + h,
           colour = "black", linewidth = 0.4)
}

gl <- function(dd, yv, lw) geom_line(data = dd,
        aes(gx, .data[[yv]], colour = trt, group = interaction(trt, chr)),
        linewidth = lw)

split_panel <- function(d, yv, brk, lo, hi, lo_brk, hi_brk, ylab, tag,
                        lw = 0.28, xaxis = FALSE) {
  if (is.na(brk)) {          # unbroken: one panel over the whole range
    p <- ggplot() + deco(xaxis) + gl(d, yv, lw) +
      scale_colour_manual(values = COL, drop = FALSE) +
      scale_y_continuous(breaks = lo_brk, expand = expansion(c(0, 0.03))) +
      coord_cartesian(xlim = c(0, xmax_all), ylim = lo) +
      labs(y = ylab, tag = tag) + theme(plot.margin = margin(7, 4, 2, 16))
    return(list(top = NULL, bot = p))
  }
  top <- ggplot() + deco() + gl(mask(d, yv, TRUE, brk), yv, lw) +
    scale_colour_manual(values = COL, drop = FALSE) +
    scale_y_continuous(breaks = hi_brk, expand = expansion(0)) +
    coord_cartesian(xlim = c(0, xmax_all), ylim = hi, clip = "off") + slash(hi, hi[1], 1) +
    labs(y = NULL, tag = tag) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), plot.margin = margin(7, 4, 2, 16))
  bot <- ggplot() + deco(xaxis) + gl(mask(d, yv, FALSE, brk), yv, lw) +
    scale_colour_manual(values = COL, drop = FALSE) +
    scale_y_continuous(breaks = lo_brk, expand = expansion(0)) +
    coord_cartesian(xlim = c(0, xmax_all), ylim = lo, clip = "off") + slash(lo, lo[2], 3) +
    labs(y = ylab) + theme(plot.margin = margin(2, 4, 2, 16))
  list(top = top, bot = bot)
}

# One point per series, as in Figure 1. geom_line warns that each group has a
# single observation and draws nothing, which is the point -- only the legend is
# wanted. Two points per group would draw a visible line across the strip.
key <- function(lev) {
  ggplot(tibble(x = 1, y = 1, g = factor(lev, levels = lev)), aes(x, y, colour = g)) +
    geom_line() + scale_colour_manual(values = COL[lev], name = NULL, drop = FALSE) +
    guides(colour = guide_legend(nrow = 1, keywidth = unit(9, "pt"))) +
    theme_void(base_size = 8) +
    theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
          legend.text = element_text(size = 7))
}

# ── build one figure from a subset of the series ─────────────────────────────
build <- function(d, lev, out) {
  d <- d %>% filter(trt %in% lev) %>% mutate(trt = factor(as.character(trt), levels = lev))
  pA <- split_panel(d, "wald", BRK_A, LO_A, HI_A, TICK_LO_A, TICK_HI_A,
                    expression(-log[10] * italic(P)), "A")
  pB <- split_panel(d, "h2", BRK_B, LO_B, HI_B, TICK_LO_B, TICK_HI_B,
                    expression(italic(h)^2), "B", xaxis = TRUE)

  parts <- list(); hts <- numeric(0)
  add <- function(p, h) { if (!is.null(p)) { parts[[length(parts)+1]] <<- p; hts <<- c(hts, h) } }
  add(pA$top, 1); add(pA$bot, 3)
  add(pB$top, 1); add(pB$bot, 3)
  add(key(lev), 0.45)

  p <- wrap_plots(parts, ncol = 1, heights = hts)
  png(out, width = W_IN, height = H_IN, units = "in", res = DPI)
  # the key's single-observation warning and the masked-out-of-range notes are
  # both expected; they would otherwise bury anything real
  suppressWarnings(print(p)); dev.off()
  cat("wrote", out, "--", length(lev), "series\n")
}

build(scans, levels(scans$trt)[levels(scans$trt) %in% unique(as.character(scans$trt))], OUT_A)
build(scans, intersect(FEMALES, unique(as.character(scans$trt))), OUT_B)
