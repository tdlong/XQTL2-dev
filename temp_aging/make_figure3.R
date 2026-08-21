# make_figure3.R — Figure 3: the SNP scan, one panel per treatment.
#
#   Rscript temp_aging/make_figure3.R
#
# DATA: the 10-replicate AGE_SY dataset -- replicates 8 and 9 dropped, those
# being the May 2023 cage (helpfiles/AGE_2024/population_assignment.txt).
# Needs process/AGE_SY/AGE_SY_4snpscan_no89.txt.gz, written on HPC3 by
# scripts_oneoffs/AGE_SY/nov_only/gather_snp.R and fetched with scp.
#
# Four panels, one per treatment, sharing Figure 1's conventions: chromosome arms
# concatenated left to right, dotted lines between arms, grey bands over
# heterochromatin, one horizontal axis across all four.
#
# The SNP scan tests every SNP, so each panel is a dense scatter rather than a
# line. It is not independent of the haplotype scan in Figure 1 -- both come from
# the same smoothed haplotype estimates.
#
# Broken y axis, built the same way as Figure 1 panel A and for the same reason:
# the chr3L locus alone sets the scale otherwise. Each panel is two ggplots
# stacked 3:1. Two details matter and are easy to lose --
#   expand = expansion(0)  the two sub-panel ranges must be mutually exclusive
#     and exhaustive, or the default 5% padding overlaps them and points in the
#     overlap are drawn TWICE, once in each.
#   mask()                 the break marks need clip = "off", which also stops
#     out-of-range points being clipped, so they spill across the figure. Setting
#     the wrong side to NA keeps them from being drawn at all.

suppressMessages({library(tidyverse); library(patchwork)})

SNP  <- "process/AGE_SY/AGE_SY_4snpscan_no89.txt.gz"
REMOTE <- "tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev"
OUT  <- "temp_aging/Figure3_plot.png"
W_IN <- 7.5; H_IN <- 6.5; DPI <- 300
PT_SIZE <- 0.35; PT_ALPHA <- 0.6

CHRS   <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
CHRLAB <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")
# dm6 euchromatin boundaries (Huynh et al. 2023 PLoS Genet 19:e1010439, Table S2)
HET <- tribble(~chr, ~eu_start, ~eu_end,
  "chrX", 2.5, 21.2, "chr2L", 0.5, 22.9, "chr2R", 1.3, 25.1,
  "chr3L", 0.7, 24.0, "chr3R", 4.5, 32.0)

# Figure 1's colours
TRT_LEV <- c("SY10 female", "SY20 female", "SY10 male", "SY20 male")
TRT_COL <- c("SY10 female" = "#F49AC2", "SY20 female" = "#D62728",
             "SY10 male"   = "#8EC7E8", "SY20 male"   = "#1F4E9C")
# Two panels, both diets overlaid within a sex. The two diets track each other
# closely, so overlaying them shows where they diverge instead of asking the
# reader to compare across panels.
PANELS <- list(A = c("SY10 female", "SY20 female"),
               B = c("SY10 male",   "SY20 male"))
PANEL_TITLE <- c(A = "females", B = "males")
# the paler diet drawn first so the darker is not buried under it
DRAW_ORDER <- c("SY10 female", "SY20 female", "SY10 male", "SY20 male")

# Fetch it if it is not here. run_snp_scans.sh chains the gather on the cluster,
# so by the time those jobs are done the file exists there; there is no reason to
# make fetching it a separate thing to remember.
if (!file.exists(SNP)) {
  cat("not found locally, fetching:\n  ", SNP, "\n", sep = "")
  dir.create(dirname(SNP), showWarnings = FALSE, recursive = TRUE)
  rc <- system2("scp", c(file.path(REMOTE, SNP), dirname(SNP)))
  if (rc != 0 || !file.exists(SNP))
    stop("scp failed. Either the gather has not finished on HPC3, or fetch it by hand:\n",
         "  scp ", file.path(REMOTE, SNP), " ", dirname(SNP), "/", call. = FALSE)
  cat("fetched", format(structure(file.size(SNP), class = "object_size"),
                        units = "auto"), "\n\n")
}

d <- read.table(SNP, header = TRUE, sep = "\t") %>% as_tibble() %>%
  filter(chr %in% CHRS) %>%
  mutate(chr = factor(chr, levels = CHRS),
         trt = factor(paste0(sugar, " ", ifelse(sex == "F", "female", "male")),
                      levels = TRT_LEV))

cat("SNPs per treatment, and the tallest:\n\n")
d %>% group_by(trt) %>%
  summarise(snps = n(), max_wald = round(max(Wald_log10p, na.rm = TRUE), 1),
            `% > 5` = round(100 * mean(Wald_log10p > 5, na.rm = TRUE), 1),
            .groups = "drop") %>% as.data.frame() %>% print(row.names = FALSE)
cat("\n")

# ── one concatenated x-axis shared by all four panels ────────────────────────
lens <- d %>% group_by(chr) %>% summarise(len = max(pos), .groups = "drop") %>%
  arrange(chr) %>% mutate(offset = lag(cumsum(len), default = 0))
d <- d %>% left_join(lens %>% select(chr, offset), by = "chr") %>%
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

# One scale across all four, split at 40. Below the break is where every panel
# has its bulk; above it is chr3L and little else. Shared rather than free so the
# four stay comparable.
BRK <- 40
LO  <- c(0, BRK);  LO_BRK <- c(0, 10, 20, 30)
HI  <- c(BRK, 170); HI_BRK <- c(50, 100, 150)

mask <- function(x, keep_above) { x[if (keep_above) x < BRK else x > BRK] <- NA_real_; x }

# one stroke marking where a sub-panel's y axis terminates; rel is the
# sub-panel's relative height, so the stroke is the same size on the page in
# both and therefore the same angle
slash <- function(yr, at, rel) {
  h <- diff(yr) * 0.030 / rel
  w <- xmax_all * 0.004
  annotate("segment", x = -w, xend = w, y = at - h, yend = at + h,
           colour = "black", linewidth = 0.4)
}

deco <- function() list(
  geom_rect(data = het_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey88", colour = NA, inherit.aes = FALSE),
  geom_vline(xintercept = chr_edges, linetype = "dotted",
             colour = "grey35", linewidth = 0.3),
  scale_x_continuous(expand = expansion(0), breaks = chr_breaks,
                     labels = CHRLAB[CHRS]),
  theme_classic(base_size = 8),
  theme(axis.title.x = element_blank(),
        plot.tag = element_text(size = 10, face = "bold"),
        plot.tag.position = c(0.004, 0.88),
        legend.position = "none"))

pts <- function(lvs, keep_above) {
  lapply(intersect(DRAW_ORDER, lvs), function(lv) {
    dd <- d %>% filter(trt == lv) %>% mutate(y = mask(Wald_log10p, keep_above))
    geom_point(data = dd, aes(gx, y), colour = TRT_COL[[lv]],
               size = PT_SIZE, alpha = PT_ALPHA, stroke = 0)
  })
}

key <- function(lvs) list(
  geom_point(data = tibble(gx = -Inf, y = -Inf, g = factor(lvs, levels = lvs)),
             aes(gx, y, colour = g), size = 1.4, na.rm = TRUE),
  scale_colour_manual(values = TRT_COL[lvs], name = NULL, drop = FALSE),
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))),
  theme(legend.position = "inside", legend.position.inside = c(0.075, 0.72),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.height = unit(8, "pt"), legend.text = element_text(size = 6.5),
        legend.margin = margin(1, 2, 1, 2)))

panel <- function(lv, tag, xaxis) {
  top <- ggplot() + deco() + pts(lv, TRUE) +
    scale_y_continuous(breaks = HI_BRK, expand = expansion(0)) +
    coord_cartesian(xlim = c(0, xmax_all), ylim = HI, clip = "off") +
    slash(HI, HI[1], 1) +
    labs(y = NULL, tag = tag) +
    annotate("text", x = xmax_all * 0.5, y = Inf, label = PANEL_TITLE[[tag]],
             vjust = 1.5, size = 2.8, fontface = "bold") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), plot.margin = margin(6, 4, 2, 16))
  bot <- ggplot() + deco() + pts(lv, FALSE) + key(lv) +
    scale_y_continuous(breaks = LO_BRK, expand = expansion(0)) +
    coord_cartesian(xlim = c(0, xmax_all), ylim = LO, clip = "off") +
    slash(LO, LO[2], 3) +
    labs(y = expression(-log[10] * italic(P))) +
    theme(plot.margin = margin(2, 4, 2, 16)) +
    (if (xaxis) NULL else theme(axis.text.x = element_blank(),
                                axis.ticks.x = element_blank()))
  list(top, bot)
}

tags <- names(PANELS)
cells <- pmap(list(PANELS, tags, tags == tags[length(tags)]), panel)
ps <- unlist(cells, recursive = FALSE)

png(OUT, width = W_IN, height = 4.4, units = "in", res = DPI)
suppressWarnings(print(wrap_plots(ps, ncol = 1,
                                  heights = rep(c(1, 3), length(PANELS)))))
dev.off()
cat("wrote", OUT, "\n")
