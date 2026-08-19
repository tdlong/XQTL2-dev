# make_figureY.R — Wald profiles of the two source cages, one panel per food.
#
# Three stacked panels sharing one concatenated x-axis and one y scale:
#   A  SY10, females
#   B  SY20, females
#   C  lab food (the Aug 2024 pilot), females
#
# Two lines per panel: the May 2023 cage and the November 2023 cage. Both were
# built from the same eight founders and maintained separately
# (helpfiles/AGE_2024/population_assignment.txt).
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript temp_aging/make_figureY.R
#
# The replicate counts are not matched -- SY10 and SY20 are May 2 against Nov 10,
# lab is May 4 against Nov 2 -- and the Wald statistic grows with them, so the raw
# profiles are not comparable in height. Each is therefore divided by its own
# replicate count, which puts all six in the same 3-6 range and downweights the
# well-replicated scans. The counts stay in the legends regardless.
#
# The y axis is deliberately low: the peaks run off the top so the low-level
# structure across the rest of the genome is legible. Raise YMAX to see them.

suppressMessages({library(tidyverse); library(patchwork)})

GATHER <- "process/AGE_2024/AGE_2024_vs_AGE_SY.txt.gz"
OUT    <- "temp_aging/FigureY_plot.png"
W_IN <- 7.5; H_IN <- 5.4; DPI <- 300
NORMALIZE <- TRUE          # divide each scan's Wald by its replicate count
YMAX <- 3.5                # ceiling; peaks above it are clipped, not dropped
YBRK <- c(0, 1, 2, 3)

CHRS   <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
CHRLAB <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")
# dm6 euchromatin boundaries (Huynh et al. 2023 PLoS Genet 19:e1010439, Table S2)
HET <- tribble(~chr, ~eu_start, ~eu_end,
  "chrX", 2.5, 21.2, "chr2L", 0.5, 22.9, "chr2R", 1.3, 25.1,
  "chr3L", 0.7, 24.0, "chr3R", 4.5, 32.0)

# cage colours, deliberately unlike the food/sex palette of Figures 1 and X
POP_COL <- c(May = "#1B7837", Nov = "#762A83")
# panel order, and how many replicates each scan actually used
PANELS <- tribble(~diet,  ~tag, ~title,
                  "SY10", "A",  "SY10",
                  "SY20", "B",  "SY20",
                  "lab",  "C",  "lab food")
NREP <- tribble(~diet, ~pop, ~n,
  "SY10","May",2, "SY10","Nov",10,
  "SY20","May",2, "SY20","Nov",10,
  "lab", "May",4, "lab", "Nov",2)

# ── data ─────────────────────────────────────────────────────────────────────
if (!file.exists(GATHER)) stop("missing ", GATHER)
g <- read.table(GATHER, header = TRUE, sep = "\t") %>% as_tibble()

d <- g %>% filter(reps %in% c("May", "Nov"), chr %in% CHRS) %>%
  transmute(chr = factor(chr, levels = CHRS), pos, diet,
            pop = factor(reps, levels = c("May", "Nov")), wald = Wald_log10p) %>%
  left_join(NREP, by = c("diet", "pop" = "pop")) %>%
  mutate(wald = if (NORMALIZE) wald / n else wald)

missing <- setdiff(paste(NREP$diet, NREP$pop), paste(d$diet, d$pop) %>% unique())
if (length(missing)) cat("NOT FOUND:", paste(missing, collapse = ", "), "\n\n")

cat(if (NORMALIZE) "Wald / n_replicates" else "Wald", " by food and cage:\n\n", sep = "")
d %>% group_by(diet, pop, n) %>%
  summarise(max = round(max(wald), 2), median = round(median(wald), 2),
            `% clipped` = round(100 * mean(wald > YMAX), 1), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n")

# ── one concatenated x-axis shared by all three panels ───────────────────────
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

# ── panels ───────────────────────────────────────────────────────────────────
panel <- function(dt, tag, title, xaxis) {
  dd <- d %>% filter(diet == dt)
  # replicate count into the key, so it travels with the line
  lab <- NREP %>% filter(diet == dt) %>%
    mutate(txt = sprintf("%s 2023  (n = %d)", pop, n))
  key <- setNames(lab$txt, lab$pop)
  dd <- dd %>% mutate(lg = factor(key[as.character(pop)], levels = key))
  cols <- setNames(POP_COL[names(key)], key)

  p <- ggplot(dd) +
    geom_rect(data = het_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey88", colour = NA, inherit.aes = FALSE) +
    geom_vline(xintercept = chr_edges, linetype = "dotted",
               colour = "grey35", linewidth = 0.3) +
    geom_line(aes(gx, wald, colour = lg, group = interaction(lg, chr)),
              linewidth = 0.3) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_x_continuous(expand = expansion(0), breaks = chr_breaks,
                       labels = CHRLAB[CHRS]) +
    scale_y_continuous(breaks = YBRK, expand = expansion(c(0, 0.04))) +
    coord_cartesian(xlim = c(0, xmax_all), ylim = c(0, YMAX)) +
    labs(y = if (NORMALIZE) expression(-log[10] * italic(P) / n)
             else expression(-log[10] * italic(P)), tag = tag) +
    annotate("text", x = xmax_all * 0.5, y = Inf, label = title,
             vjust = 1.6, size = 2.8, fontface = "bold") +
    theme_classic(base_size = 8) +
    theme(axis.title.x = element_blank(),
          plot.tag = element_text(size = 10, face = "bold"),
          plot.tag.position = c(0.004, 0.92),
          legend.position = "inside",
          legend.position.inside = c(0.085, 0.80),
          legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
          legend.key.width = unit(11, "pt"), legend.key.height = unit(8, "pt"),
          legend.text = element_text(size = 6.5),
          plot.margin = margin(5, 4, 2, 16))
  if (!xaxis) p <- p + theme(axis.text.x = element_blank(),
                             axis.ticks.x = element_blank())
  p
}

ps <- pmap(list(PANELS$diet, PANELS$tag, PANELS$title,
                PANELS$diet == PANELS$diet[nrow(PANELS)]), panel)

png(OUT, width = W_IN, height = H_IN, units = "in", res = DPI)
suppressWarnings(print(wrap_plots(ps, ncol = 1)))
dev.off()
cat("wrote", OUT, "\n")
