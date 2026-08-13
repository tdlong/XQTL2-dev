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

suppressMessages({library(tidyverse); library(ggbreak); library(aplot)})

LONG      <- "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
VARCOMP   <- "process/AGE_SY_splithalf/H2_varcomp_by_window.txt.gz"
LIVE_SCAN_DIR <- NULL          # e.g. "process/AGE_SY/Scans" for the 12-rep scans
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

scans <- if (is.null(LIVE_SCAN_DIR)) {
  long %>% group_by(chr, pos, sugar, sex) %>%
    summarise(wald = mean(Wald_log10p), h2 = mean(Cutl_H2), .groups = "drop")
} else {
  expand_grid(sugar = c("SY10","SY20"), sex = c("F","M")) %>%
    mutate(scan = paste0("AGE_", sugar, "_", sex)) %>%
    pmap_dfr(function(sugar, sex, scan) {
      read.table(file.path(LIVE_SCAN_DIR, scan, paste0(scan, ".scan.txt")),
                 header = TRUE) %>%
        transmute(chr, pos = as.integer(pos), sugar, sex,
                  wald = Wald_log10p, h2 = Cutl_H2)
    })
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
# Broken y axis on every panel via ggbreak, so the chr3L locus does not set the
# scale for the whole genome. scales = 1/3 makes the upper segment a third the
# height of the lower, i.e. the 3:1 split.

deco <- function(tag, xaxis = FALSE) {
  th <- theme_classic(base_size = 8) +
    theme(axis.title.x = element_blank(),
          plot.tag = element_text(size = 10, face = "bold"),
          plot.tag.position = c(0.004, 0.97),
          plot.margin = margin(2, 4, 1, 2),
          legend.position = "none",
          # ggbreak mirrors the axis on the right; drop it
          axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(),
          axis.line.y.right = element_blank(), axis.title.y.right = element_blank())
  if (!xaxis) th <- th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  list(
    geom_rect(data = het_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey88", colour = NA, inherit.aes = FALSE),
    geom_vline(xintercept = chr_edges, linetype = "dotted", colour = "grey35", linewidth = 0.3),
    scale_x_continuous(limits = c(0, xmax_all), expand = expansion(0),
                       breaks = chr_breaks, labels = CHRLAB[CHRS]),
    labs(tag = tag), th)
}

pA <- ggplot(scans, aes(gx, wald, colour = trt)) +
  geom_line(linewidth = 0.28) +
  scale_colour_manual(values = TRT_COL) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 25, 50, 75, 100)) +
  scale_y_break(c(15, 17), scales = 1/3, space = 0.12, ticklabels = c(25, 50, 75, 100)) +
  labs(y = expression(-log[10] * italic(P))) +
  deco("A")

pB <- ggplot(scans, aes(gx, h2, colour = trt)) +
  geom_line(linewidth = 0.28) +
  scale_colour_manual(values = TRT_COL) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 3, 4, 5)) +
  scale_y_break(c(2, 2.2), scales = 1/3, space = 0.12, ticklabels = c(3, 4, 5)) +
  labs(y = expression(italic(h)^2)) +
  deco("B")

pC <- ggplot(cmp, aes(gx, y, colour = term)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
  geom_line(linewidth = 0.35) +
  scale_colour_manual(values = CMP_COL) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.25, 0.5, 0.75, 1, 2)) +
  scale_y_break(c(0.75, 0.82), scales = 1/3, space = 0.12, ticklabels = c(1, 2)) +
  labs(y = expression(italic(h)^2)) +
  deco("C", xaxis = TRUE)

key <- function(lev, col) {
  ggplot(tibble(x = 1, y = 1, g = factor(lev, levels = lev)), aes(x, y, colour = g)) +
    geom_line() + scale_colour_manual(values = col, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(nrow = 1, keywidth = unit(9, "pt"),
                                 override.aes = list(linewidth = 0.8))) +
    theme_void(base_size = 8) + theme(legend.position = "bottom")
}
legs <- ggplotify::as.ggplot(
  cowplot::plot_grid(cowplot::get_legend(key(TRT_LEV, TRT_COL)),
                     cowplot::get_legend(key(CMP_LEV, CMP_COL)),
                     nrow = 1, rel_widths = c(2.1, 1)))

g <- aplot::plot_list(pA, pB, pC, legs, ncol = 1, heights = c(1, 1, 1.12, 0.1))

png(OUT, width = W_IN, height = H_IN, units = "in", res = DPI)
print(g)
invisible(dev.off())
cat("Wrote:", OUT, "\n")
