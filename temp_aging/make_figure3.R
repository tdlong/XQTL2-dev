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

suppressMessages({library(tidyverse); library(patchwork)})

SNP  <- "process/AGE_SY/AGE_SY_4snpscan_no89.txt.gz"
OUT  <- "temp_aging/Figure3_plot.png"
W_IN <- 7.5; H_IN <- 6.5; DPI <- 300
PT_SIZE <- 0.05; PT_ALPHA <- 0.25

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
TAGS <- c("A", "B", "C", "D")

if (!file.exists(SNP)) stop("missing ", SNP,
  "\nRun scripts_oneoffs/AGE_SY/nov_only/gather_snp.R on HPC3 and scp it back.",
  call. = FALSE)

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

# one y scale across the four, so panels are comparable
YMAX <- max(d$Wald_log10p, na.rm = TRUE) * 1.02

panel <- function(lv, tag, xaxis) {
  ggplot(d %>% filter(trt == lv)) +
    geom_rect(data = het_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey88", colour = NA, inherit.aes = FALSE) +
    geom_vline(xintercept = chr_edges, linetype = "dotted",
               colour = "grey35", linewidth = 0.3) +
    geom_point(aes(gx, Wald_log10p), colour = TRT_COL[[lv]],
               size = PT_SIZE, alpha = PT_ALPHA, stroke = 0) +
    scale_x_continuous(expand = expansion(0), breaks = chr_breaks,
                       labels = CHRLAB[CHRS]) +
    scale_y_continuous(expand = expansion(c(0, 0.04))) +
    coord_cartesian(xlim = c(0, xmax_all), ylim = c(0, YMAX)) +
    labs(y = expression(-log[10] * italic(P)), tag = tag) +
    annotate("text", x = xmax_all * 0.5, y = Inf, label = as.character(lv),
             vjust = 1.6, size = 2.8, fontface = "bold") +
    theme_classic(base_size = 8) +
    theme(axis.title.x = element_blank(),
          plot.tag = element_text(size = 10, face = "bold"),
          plot.tag.position = c(0.004, 0.92),
          legend.position = "none",
          plot.margin = margin(5, 4, 2, 16)) +
    (if (xaxis) NULL else theme(axis.text.x = element_blank(),
                                axis.ticks.x = element_blank()))
}

ps <- pmap(list(TRT_LEV, TAGS, TRT_LEV == TRT_LEV[length(TRT_LEV)]), panel)

png(OUT, width = W_IN, height = H_IN, units = "in", res = DPI)
suppressWarnings(print(wrap_plots(ps, ncol = 1)))
dev.off()
cat("wrote", OUT, "\n")
