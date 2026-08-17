# make_figure2.R — seven zoomed peaks, plus a legend/control cell.
#
#   Rscript temp_aging/make_figure2.R
#
# Needs process/AGE_SY/AGE_SY_zoom_means.txt.gz (make_zoom_means.R, run on HPC3)
# and process/AGE_SY/AGE_SY_4scan.txt.gz.
#
# LAYOUT. Four columns by two rows of cells. Seven cells are a peak; the eighth
# is a legend over the control panel. Each peak cell is two sub-panels:
#
#   upper 1/3   -log10 P for all four treatments, over 1 Mb centred on the peak
#   lower 2/3   change in founder frequency (Z - C) for the treatment with the
#               highest Wald at that peak, founders in fixed colours, founders
#               averaging under 2.5% in the controls faded (the XQTL2.Xplore
#               convention, github.com/tdlong/XQTL2.Xplore, XQTL_change_average)
#
# THE CONTROL PANEL asks whether these hits are just founders being purged from
# the cages over time. At each peak we take the founder that drops most under
# selection -- excluding the faded rare ones, whose changes are noise -- and plot
# its frequency in the CONTROL pools against replicate, replicates having accrued
# over months. Seven lines, one per peak. Lines that are flat say the founder is
# not on its way out of the cage and the drop is selection, not drift or purging.

suppressMessages({library(tidyverse); library(patchwork)})

MEANS <- "process/AGE_SY/AGE_SY_zoom_means.txt.gz"
SCAN  <- "process/AGE_SY/AGE_SY_4scan.txt.gz"
OUT   <- "temp_aging/Figure2_plot.png"
WIN   <- 0.5e6                      # half-width actually plotted
RARE  <- 0.025                      # founders below this in controls are faded

# Common y scales, so panels are comparable at a glance. chr3L is the exception
# in both sub-panels -- it is an order of magnitude bigger and would flatten
# everything else. Limits sit just outside the breaks because two panels would
# otherwise clip: chr2L:10.71 reaches Wald 30.5 and chr2L:11.95 reaches 0.061.
BIG   <- "chr3L:9.31"
YLIM  <- list(
  wald = list(def = list(lim = c(0, 31),           brk = c(0, 10, 20, 30)),
              big = list(lim = c(0, 212),          brk = c(0, 50, 100, 150, 200))),
  freq = list(def = list(lim = c(-0.062, 0.062),   brk = c(-0.06, -0.03, 0, 0.03, 0.06)),
              # chr3L has its own scale anyway, so it is trimmed to its data:
              # the trough reaches -0.110 and nothing rises past 0.046.
              big = list(lim = c(-0.115, 0.052),   brk = c(-0.10, -0.05, 0, 0.05))))
yspec <- function(which, lc) YLIM[[which]][[if (lc == BIG) "big" else "def"]]

if (!file.exists(MEANS)) stop("missing ", MEANS,
  "\nRun temp_aging/make_zoom_means.R on HPC3 and scp the result back.")

TRT_LEV <- c("SY10 female", "SY20 female", "SY10 male", "SY20 male")
TRT_COL <- c("SY10 female" = "#F49AC2", "SY20 female" = "#D62728",
             "SY10 male"   = "#8EC7E8", "SY20 male"   = "#1F4E9C")
trt_lab <- function(sugar, sex) factor(paste0(sugar, " ", ifelse(sex == "F", "female", "male")),
                                       levels = TRT_LEV)

# Founder colours: the XQTL2.Xplore palette, fixed so a founder is the same
# colour in every panel.
FOUNDER_COL <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#D55E00", "#CC79A7", "#000000", "#990099")

# Peak colours. The control panel has one line per peak, not per founder -- two
# peaks can share a top-dropping founder. Each cell's locus label is drawn in its
# peak colour, which is the key, so the control panel needs no legend of its own.
PEAK_COL <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
              "#66A61E", "#E6AB02", "#A6761D")

means <- read.table(MEANS, header = TRUE, sep = "\t") %>% as_tibble()
scan  <- read.table(SCAN,  header = TRUE, sep = "\t") %>% as_tibble()

founders <- sort(unique(means$founder))
names(FOUNDER_COL) <- founders

loci <- means %>% distinct(locus, chr, peak_pos) %>%
  arrange(match(chr, c("chrX","chr2L","chr2R","chr3L")), peak_pos)
names(PEAK_COL) <- loci$locus

# treatment with the highest Wald at each peak -- the one whose frequency change
# gets drawn
best <- loci %>% pmap_dfr(function(locus, chr, peak_pos)
  scan %>% filter(chr == !!chr, pos == !!peak_pos) %>% slice_max(Wald_log10p, n = 1) %>%
    transmute(locus, sugar, sex, wald = Wald_log10p))

# ---- per-peak frequency change, averaged over replicates ------------------
dfreq <- means %>%
  inner_join(best %>% select(locus, sugar, sex), by = c("locus","sugar","sex")) %>%
  pivot_wider(names_from = TRT, values_from = freq, names_prefix = "f") %>%
  group_by(locus, chr, peak_pos, pos, founder) %>%
  summarise(Dfreq = mean(fZ - fC, na.rm = TRUE),
            fC    = mean(fC, na.rm = TRUE), .groups = "drop") %>%
  group_by(locus, founder) %>% mutate(rare = mean(fC, na.rm = TRUE) < RARE) %>%
  ungroup() %>% mutate(alpha = ifelse(rare, 0.2, 1)) %>%
  arrange(alpha, founder)                      # faded lines drawn behind

# ---- the founders that move most at each peak, ignoring the faded ones -----
# One that rises under selection and one that falls, per peak.
at_peak <- dfreq %>% filter(!rare, pos == peak_pos)
movers <- bind_rows(
  at_peak %>% group_by(locus) %>% slice_max(Dfreq, n = 1) %>% mutate(dir = "up"),
  at_peak %>% group_by(locus) %>% slice_min(Dfreq, n = 1) %>% mutate(dir = "down")) %>%
  ungroup() %>% select(locus, founder, dir, Dfreq_at_peak = Dfreq)

# their frequency in the controls, by replicate
ctrl <- means %>%
  inner_join(best %>% select(locus, sugar, sex), by = c("locus","sugar","sex")) %>%
  inner_join(movers %>% select(locus, founder, dir), by = c("locus","founder")) %>%
  filter(TRT == "C", pos == peak_pos) %>%
  transmute(locus, founder, dir, REP, freq)

# ---- panels ---------------------------------------------------------------
base <- theme_bw(7) + theme(panel.grid.minor = element_blank(),
          plot.margin = margin(1, 2, 1, 2), legend.position = "none",
          plot.title = element_blank())

# The control panels reserve no space for a sign. An earlier version padded them
# to the width of "-0.06" so the two boxes would match, which was only needed
# while cells 4 and 8 were nested and sized independently; in one flat layout
# patchwork aligns the axis block across the whole column, so one decimal and no
# padding is enough.
ylab_zoom <- function(x) sprintf("%.2f", x)

# Only panels that introduce a scale carry y tick labels. Repeating "-0.06" on
# every frequency panel is not just clutter: patchwork sizes a column's axis
# block to its widest label, so cell 4's labels were forcing two characters of
# empty buffer to the left of the control panels' own "0.6".
noy <- function(show) if (show) NULL else
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ylab_ctrl <- function(x) sprintf("%.1f", x)

wald_panel <- function(lc, ylab = NULL, showy = FALSE) {
  L <- loci %>% filter(locus == lc)
  d <- scan %>% filter(chr == L$chr, abs(pos - L$peak_pos) <= WIN) %>%
       mutate(trt = trt_lab(sugar, sex), mb = pos/1e6)
  ggplot(d, aes(mb, Wald_log10p, colour = trt)) +
    geom_line(linewidth = 0.35) +
    scale_colour_manual(values = TRT_COL) +
    scale_x_continuous(expand = expansion(0)) +
    scale_y_continuous(limits = yspec("wald", lc)$lim,
                       breaks = yspec("wald", lc)$brk, expand = expansion(0)) +
    labs(x = NULL, y = ylab) + base + noy(showy) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    annotate("text", x = -Inf, y = Inf, label = lc, hjust = -0.08, vjust = 1.4,
             size = 2.1, colour = PEAK_COL[[lc]], fontface = "bold")
}

freq_panel <- function(lc, ylab = NULL, showy = FALSE) {
  d <- dfreq %>% filter(locus == lc) %>% mutate(mb = pos/1e6)
  ggplot(d, aes(mb, Dfreq, colour = founder, alpha = alpha,
                group = interaction(founder, locus))) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.2) +
    geom_line(linewidth = 0.4) +
    scale_colour_manual(values = FOUNDER_COL) + scale_alpha_identity() +
    scale_x_continuous(expand = expansion(0)) +
    scale_y_continuous(limits = yspec("freq", lc)$lim,
                       breaks = yspec("freq", lc)$brk, expand = expansion(0),
                       labels = ylab_zoom) +
    labs(x = NULL, y = ylab) + base + noy(showy)
}

# Founder key, drawn as a plot rather than an extracted legend: cowplot's
# get_legend returns a gtable, which patchwork will not compose. The treatment
# colours are not repeated here -- they are defined in Figure 1 -- and the peak
# colours need no key, each cell's locus label being drawn in its own colour.
founder_legend <- function() {
  fk <- tibble(lab = founders, x = seq_along(founders))
  ggplot(fk) +
    geom_point(aes(x, 1, colour = lab), size = 1.6) +
    geom_text(aes(x + 0.09, 1, label = lab), hjust = 0, size = 2.1) +
    scale_colour_manual(values = FOUNDER_COL) +
    scale_x_continuous(limits = c(0.6, length(founders) + 0.8)) +
    theme_void(7) + theme(legend.position = "none",
                          plot.margin = margin(2, 3, 0, 3))
}

# The two halves of the eighth cell. At each peak, the founder enriched most in
# the long-lived pool (most protective) and the one depleted most (most
# susceptible), faded rare founders excluded, shown as their frequency in the
# CONTROLS across replicates. If these hits were just founders on their way out
# of the cage, the susceptible set would slope down here too.
control_panel <- function(d, lab, show_x) {
  ggplot(ctrl %>% filter(dir == d), aes(REP, freq, colour = locus, group = locus)) +
    geom_line(linewidth = 0.4) + geom_point(size = 0.5) +
    scale_colour_manual(values = PEAK_COL) +
    scale_x_continuous(breaks = c(1, 4, 8, 12), expand = expansion(0.03)) +
    scale_y_continuous(limits = c(0, 0.75), breaks = c(0, 0.2, 0.4, 0.6),
                       expand = expansion(0), labels = ylab_ctrl) +
    labs(x = if (show_x) "replicate" else NULL, y = NULL) + base +
    (if (show_x) NULL else theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank())) +
    annotate("text", x = -Inf, y = Inf, label = lab,
             hjust = -0.06, vjust = 1.5, size = 2.1)
}

# No nested blocks. patchwork equalises panel widths only for plots that are
# direct children of one layout, so everything sits in a single five-column
# design: three zoom columns, a narrow column carrying nothing but the rotated
# title for the eighth cell, and a fourth zoom column. YGAP is nearly all the
# slack there is -- at 0.10 the gap before column 4 measures 181 px against 93
# elsewhere, and at 0.03 it is 152 px. The floor is set by the y tick labels of
# that column, not by the title. Nesting cells 4 and 8 in
# their own sub-layouts is what made their boxes different widths, because each
# then got whatever space was left after its own axis material.
YGAP  <- 0.07
ylab_el <- wrap_elements(grid::textGrob("haplotype freq in controls",
             rot = 90, gp = grid::gpar(fontsize = 7)))

lc <- loci$locus
# y titles only on the leftmost column of cells (1 and 5)
YW <- expression(-log[10]*italic(P)); YF <- expression(Delta*" frequency")
# leftmost cell of each row, plus chr3L which is on its own scale in both rows
SHOWY <- c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE)
YTITLE <- c(list(YW), rep(list(NULL), 3), list(YW), rep(list(NULL), 2))
YTITLF <- c(list(YF), rep(list(NULL), 3), list(YF), rep(list(NULL), 2))
W  <- pmap(list(lc, YTITLE, SHOWY), wald_panel)
F_ <- pmap(list(lc, YTITLF, SHOWY), freq_panel)

# one shared x title under the three columns of zoom cells; the control cell
# keeps its own, since its x is replicate and not position
# The row's axis block reserves height for the control panel's "replicate" title,
# so this strip starts well below the tick labels and "Mb" sat 64 px under them
# against 13 px for "replicate". plot.margin does not move it -- patchwork fixes
# the row -- so the text is drawn above its own panel with clipping off.
MBV <- -1.6
xlab_strip <- ggplot() + annotate("text", x = 0, y = 1, label = "Mb", size = 2.4,
                                  vjust = MBV) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

# Six unit rows per cell row, so the Wald:frequency 1:2 split lands on whole rows
# (2 and 4) while the eighth cell divides in half (3 and 3). Letters map to plots
# alphabetically against the order the plots are added below.
design <- c("ABC#D", "ABC#D",
            "EFG#H", "EFG#H", "EFG#H", "EFG#H",
            "IJKLM", "IJKLM",
            "NOPLM",
            "NOPLQ", "NOPLQ", "NOPLQ",
            "RRR##",
            "SSSSS")
p <- W[[1]] + W[[2]] + W[[3]] + W[[4]] +
     F_[[1]] + F_[[2]] + F_[[3]] + F_[[4]] +
     W[[5]] + W[[6]] + W[[7]] +
     ylab_el +
     control_panel("up", "most protective", FALSE) +
     F_[[5]] + F_[[6]] + F_[[7]] +
     control_panel("down", "most susceptible", TRUE) +
     xlab_strip + founder_legend() +
     plot_layout(design = paste(design, collapse = "\n"),
                 widths  = c(1, 1, 1, YGAP, 1),
                 heights = c(rep(1, 12), 0.20, 0.34))

png(OUT, width = 7.5, height = 5.2, units = "in", res = 300)
print(p)
dev.off()

cat("wrote", OUT, "\n\n")
cat("treatment drawn in each frequency sub-panel:\n")
best %>% mutate(wald = round(wald, 1)) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nfounders moving most at each peak (faded rare founders excluded):\n")
movers %>% mutate(Dfreq_at_peak = round(Dfreq_at_peak, 3)) %>%
  left_join(ctrl %>% group_by(locus, dir) %>%
    summarise(ctrl_first = round(freq[which.min(REP)], 3),
              ctrl_last  = round(freq[which.max(REP)], 3),
              ctrl_slope = round(coef(lm(freq ~ REP))[2], 4), .groups = "drop"),
    by = c("locus", "dir")) %>% arrange(dir, Dfreq_at_peak) %>%
  as.data.frame() %>% print(row.names = FALSE)
