# make_talk_mixing.R -- the classic block diagram of Y = X b + e, drawn to the
# dimensions of the pool-seq problem: thousands of SNPs down, eight founders
# across. Blocks are lightly filled with what they hold (Y greys, X black/white,
# b founder colours, e faint).
#
#   Rscript temp_aging/analysis/make_talk_mixing.R   ->  Talk_Yxb.png

suppressMessages({library(tidyverse); library(patchwork)})
source("temp_aging/analysis/palette.R")
BASE <- 16; set.seed(5)
fnd <- names(FOUNDER_COL)
b   <- c(AB8 = .052, B1 = .145, B2 = .031, B3 = .151, B4 = .016, B5 = .041, B6 = .539, B7 = .024)
n <- 24                                          # rows drawn; stands for thousands
X <- matrix(rbinom(n * 8, 1, .4), n, 8); X[rowSums(X) == 0, 1] <- 1
Y <- as.numeric(X %*% b) + rnorm(n, 0, .04); Y <- pmin(pmax(Y, 0), 1)
e <- Y - as.numeric(X %*% b)

blk <- function(base, w, h, title, dims, words, ylim = c(0, h)) {
  base + coord_equal(xlim = c(0, w), ylim = ylim, expand = FALSE, clip = "off") + theme_void(BASE) +
    annotate("rect", xmin = 0, xmax = w, ymin = 0, ymax = h, fill = NA, colour = "black", linewidth = 1) +
    annotate("text", x = w/2, y = h + 2.2, label = title, size = 10, fontface = "bold") +
    annotate("text", x = w/2, y = h + .8, label = dims, size = 4.5, colour = "grey30") +
    annotate("text", x = w/2, y = -1.0, label = words, size = 4.6, lineheight = .95, vjust = 1) +
    theme(plot.margin = margin(40, 4, 80, 4))
}
Yg <- ggplot(tibble(y = n:1 - .5, v = Y), aes(1, y, fill = v)) + geom_tile(width = 1.9, height = .92) +
  scale_fill_gradient(low = "grey96", high = "grey35", limits = c(0, 1), guide = "none")
Xg <- ggplot(as_tibble(X) %>% mutate(y = n:1 - .5) %>% pivot_longer(-y, names_to = "f", values_to = "a") %>%
             mutate(x = as.integer(sub("V", "", f)) - .5), aes(x, y, fill = factor(a))) + geom_tile(width = .92, height = .92) +
  scale_fill_manual(values = c(`0` = "white", `1` = "black"), guide = "none")
bg <- ggplot(tibble(y = 8:1 - .5, f = factor(fnd, levels = fnd)), aes(1, y, fill = f)) + geom_tile(width = 1.9, height = .92) +
  geom_text(aes(label = f), size = 4, colour = ifelse(fnd %in% c("B6", "B7"), "white", "black")) +
  scale_fill_manual(values = FOUNDER_COL, guide = "none")
eg <- ggplot(tibble(y = n:1 - .5, v = e), aes(1, y, fill = v)) + geom_tile(width = 1.9, height = .92) +
  scale_fill_gradient2(low = "#1F4E9C", mid = "white", high = "#C0392B", limits = c(-.12, .12), guide = "none")

pY <- blk(Yg, 2, n, "Y", "n × 1", "SNP frequencies\nin the pool\n\nobserved;\nn in the thousands")
pX <- blk(Xg, 8, n, "X", "n × 8", "founder genotypes\nat each SNP\n\nknown; 0 or 1")
pb <- blk(bg, 2, 8, "b", "8 × 1", "founder frequencies\nin the pool\n\nunknown;\n≥ 0, sum to 1",
          ylim = c(-(n - 8)/2, n - (n - 8)/2))
pe <- blk(eg, 2, n, "e", "n × 1", "sequencing error\n\nbinomial,\n50–100 reads")
op <- function(s) ggplot() + annotate("text", 0, 0, label = s, size = 13) + theme_void() +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), clip = "off") + theme(plot.margin = margin(40, 0, 80, 0))

p <- (pY | op("=") | pX | op("×") | pb | op("+") | pe) + plot_layout(widths = c(2, 1.2, 8, 1.2, 2, 1.2, 2)) +
  plot_annotation(
    title = "Y = X b + e",
    subtitle = "Thousands of equations, eight unknowns.\nb is the least-squares solution: the founder frequencies that best reproduce the observed SNP frequencies.",
    theme = theme(plot.title = element_text(size = BASE + 8, face = "bold", hjust = .5),
                  plot.subtitle = element_text(size = BASE, hjust = .5, lineheight = 1.05),
                  plot.margin = margin(10, 20, 10, 20)))
ggsave("temp_aging/analysis/Talk_Yxb.png", p, width = 13, height = 8.5, dpi = 220, bg = "white")
cat("Wrote temp_aging/analysis/Talk_Yxb.png\n")
