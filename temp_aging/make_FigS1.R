# make_FigS1.R — Supplementary Figure 1: heritability against the Wald statistic.
#
#   Rscript temp_aging/make_FigS1.R
#
# The Methods claim that T and h2 are proportional, with the constant fixed by
# the chromosomes sampled and the selection intensity. This shows it: h2 against
# T for every window, in each of the four treatments and separately for the X,
# with the proportional relation drawn through them.
#
# Run from the repo root; needs process/AGE_SY/AGE_SY_4scan_no89.txt.gz.

suppressMessages({library(tidyverse); library(patchwork)})

SCAN <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
OUT  <- "temp_aging/FigureS1_plot.png"
stopifnot(file.exists(SCAN))

TRT_LEV <- c("SY10 female", "SY20 female", "SY10 male", "SY20 male")
TRT_COL <- c("#F48FB1", "#C62828", "#81D4FA", "#1A448E")

d <- read.table(SCAN, header = TRUE, sep = "\t") %>% as_tibble() %>%
  mutate(sex = ifelse(sex %in% c("FALSE", "F"), "F", "M"),
         trt = factor(paste0(sugar, " ", ifelse(sex == "F", "female", "male")),
                      levels = TRT_LEV),
         arm = ifelse(chr == "chrX", "X", "autosome"),
         # the chi-square behind the reported -log10 P, on m-1 = 7 df
         T   = qchisq(-Wald_log10p * log(10), df = 7, lower.tail = FALSE, log.p = TRUE)) %>%
  filter(is.finite(T), T > 0)

# one proportional constant, fitted through the origin on autosomal windows
cslope <- coef(lm(H2 ~ 0 + T, data = d %>% filter(arm == "autosome")))[[1]]
cat(sprintf("slope h2/T = %.5f  (1/slope = %.0f)\n", cslope, 1/cslope))

bin <- function(x) x %>%
  mutate(b = cut(T, breaks = c(0, 2^(seq(1, 9, by = 0.5))), labels = FALSE)) %>%
  group_by(trt, arm, b) %>%
  summarise(T = median(T), h2 = median(H2), n = n(), .groups = "drop") %>%
  filter(n >= 20)

pts <- bin(d)
line <- tibble(T = seq(0, max(pts$T) * 1.05, length.out = 100)) %>% mutate(h2 = cslope * T)

base <- function(p, xmax) p +
  geom_line(data = line %>% filter(T <= xmax), aes(T, h2), inherit.aes = FALSE,
            colour = "grey35", linewidth = 0.4, linetype = "22") +
  coord_cartesian(xlim = c(0, xmax), ylim = c(NA, cslope * xmax * 1.1)) +
  scale_colour_manual(values = setNames(TRT_COL, TRT_LEV), name = NULL) +
  labs(x = expression(Wald~statistic~italic(T)), y = expression(italic(h)^2~("%"))) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.key.width = unit(14, "pt"))

xa <- max(pts$T[pts$arm == "autosome"]) * 1.05
xx <- max(pts$T[pts$arm == "X"]) * 1.05
pA <- base(ggplot(pts %>% filter(arm == "autosome"), aes(T, h2, colour = trt)) +
             geom_point(size = 1.3), xa) +
  ggtitle("autosomes") + theme(plot.title = element_text(size = 9, face = "bold"))

pB <- base(ggplot(pts %>% filter(arm == "X"), aes(T, h2, colour = trt)) +
             geom_point(size = 1.3), xx) +
  ggtitle("X chromosome (note the axis range)") + theme(plot.title = element_text(size = 9, face = "bold"))

fig <- (pA | pB) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(OUT, fig, width = 7.0, height = 3.4, dpi = 300)
cat("wrote", OUT, "\n\n")

d %>% group_by(trt, arm) %>%
  summarise(windows = n(), `median h2/T` = signif(median(H2/T), 3), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
