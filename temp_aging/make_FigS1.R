# make_FigS1.R — Supplementary Figure 1: heritability as a function of the Wald.
#
#   Rscript temp_aging/make_FigS1.R
#
# The Methods give h2 as a function of the Wald statistic, the chromosomes
# sampled and the selection intensity. The relation is general, so this pools
# every window in every treatment: the curve is the prediction, the points are
# the h2 estimated from the data. Run from the repo root.

suppressMessages(library(tidyverse))

SCAN <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
OUT  <- "temp_aging/FigureS1_plot.png"
stopifnot(file.exists(SCAN))

d <- read.table(SCAN, header = TRUE, sep = "\t") %>% as_tibble() %>%
  mutate(T = qchisq(-Wald_log10p * log(10), df = 7, lower.tail = FALSE, log.p = TRUE)) %>%
  filter(is.finite(T), T > 0)

# h2 = c (T - df): under the null T averages df, so heritability tracks the
# non-centrality rather than T itself. One constant, fitted over all windows.
cslope <- coef(lm(H2 ~ 0 + I(T - 7), data = d))[[1]]
cat(sprintf("h2 = %.5f (T - 7)\n", cslope))

# binned on T, log-spaced: nearly every window sits at low T, so a linear axis
# stacks them all on the origin
obs <- d %>% filter(T > 7, H2 > 0) %>%
  mutate(b = cut(log10(T), breaks = seq(0.8, 3.0, by = 0.06), labels = FALSE)) %>%
  group_by(b) %>%
  summarise(T = median(T), h2 = median(H2), n = n(), .groups = "drop") %>%
  filter(n >= 50)

pred <- tibble(T = 10^seq(log10(8), log10(max(obs$T) * 1.1), length.out = 400)) %>%
  mutate(h2 = cslope * (T - 7))

# LWP is monotone in T, so the same axis carries both
lwp_ticks <- c(1, 2, 5, 7.5, 15, 40, 100)
T_at <- qchisq(-lwp_ticks * log(10), df = 7, lower.tail = FALSE, log.p = TRUE)

fig <- ggplot(obs, aes(T, h2)) +
  geom_line(data = pred, colour = "grey45", linewidth = 0.6) +
  geom_point(size = 1.5, colour = "#1A448E") +
  scale_x_log10(breaks = c(10, 20, 50, 100, 200, 500),
                sec.axis = dup_axis(breaks = T_at, labels = lwp_ticks, name = "LWP")) +
  scale_y_log10(breaks = c(0.03, 0.1, 0.3, 1, 3)) +
  annotation_logticks(sides = "bl", linewidth = 0.2,
                      short = unit(2,"pt"), mid = unit(3,"pt"), long = unit(4,"pt")) +
  labs(x = expression(Wald~statistic~italic(T)), y = expression(italic(h)^2~("%"))) +
  theme_classic(base_size = 9)

ggsave(OUT, fig, width = 4.8, height = 3.6, dpi = 300)
cat("wrote", OUT, "\n\n")
for (l in c(2, 5, 7.5, 15)) {
  Tv <- qchisq(-l * log(10), df = 7, lower.tail = FALSE, log.p = TRUE)
  cat(sprintf("  LWP %5.1f -> T %6.1f -> h2 %.2f%%\n", l, Tv, cslope * (Tv - 7)))
}
