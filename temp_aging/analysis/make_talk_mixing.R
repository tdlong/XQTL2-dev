# make_talk_mixing.R -- founder frequencies from pool-seq as a loss being
# minimised: the pie (b) changes while sum (Y - X b)^2 shrinks.
#
#   Rscript temp_aging/analysis/make_talk_mixing.R
#
#   Talk_loss.gif        the animation: gradient descent on b, projected onto
#                        the simplex each step, 60 frames
#   Talk_loss_strip.png  six frames of it side by side, for a static deck
#   Talk_loss_frames/    the individual frames
#
# Toy data: 8 founders, b from the chr3L 9.30 control frequencies, 40 SNPs,
# pool-seq noise at 80x. Started from equal founders.

suppressMessages({library(tidyverse); library(patchwork); library(magick)})
source("temp_aging/analysis/palette.R")
BASE <- 16; set.seed(5)
fnd <- names(FOUNDER_COL)
b_true <- c(AB8 = .052, B1 = .145, B2 = .031, B3 = .151, B4 = .016, B5 = .041, B6 = .539, B7 = .024)
n <- 40; COV <- 80
X <- matrix(rbinom(n * 8, 1, .45), n, 8, dimnames = list(NULL, fnd)); X[rowSums(X) == 0, 1] <- 1
Y <- rbinom(n, COV, as.numeric(X %*% b_true)) / COV
loss <- function(b) sum((Y - X %*% b)^2)

# projection onto {b >= 0, sum b = 1}
proj <- function(v) { u <- sort(v, decreasing = TRUE); css <- cumsum(u)
  rho <- max(which(u - (css - 1) / seq_along(u) > 0)); th <- (css[rho] - 1) / rho; pmax(v - th, 0) }
NIT <- 60; lr <- 0.005
path <- matrix(NA, NIT + 1, 8, dimnames = list(NULL, fnd)); path[1, ] <- rep(1/8, 8)
for (i in 1:NIT) path[i + 1, ] <- proj(path[i, ] - lr * as.numeric(-2 * t(X) %*% (Y - X %*% path[i, ])))
L <- apply(path, 1, loss)
cat(sprintf("loss: start %.3f, end %.4f; b-hat vs truth max abs error %.3f\n", L[1], L[NIT + 1], max(abs(path[NIT + 1, ] - b_true))))

frame <- function(i) {
  bi <- path[i, ]
  pie <- ggplot(tibble(f = factor(fnd, levels = fnd), b = bi), aes(1, b, fill = f)) +
    geom_col(width = 1, colour = "white", linewidth = .8) + coord_polar(theta = "y") +
    geom_text(aes(label = ifelse(b > .06, sprintf("%s", f), "")), position = position_stack(vjust = .5),
              size = 4.5, colour = ifelse(fnd %in% c("B6", "B7"), "white", "black")) +
    scale_fill_manual(values = FOUNDER_COL, guide = "none") +
    labs(title = "b: the founder frequencies", subtitle = "the unknowns -- the weights") +
    theme_void(BASE) + theme(plot.title = element_text(hjust = .5, size = BASE + 1, face = "bold"),
                             plot.subtitle = element_text(hjust = .5, size = BASE - 1))
  fit <- ggplot(tibble(s = 1:n, Y = Y, Xb = as.numeric(X %*% bi)), aes(s)) +
    geom_segment(aes(xend = s, y = Y, yend = Xb), colour = "#C0392B", linewidth = .9, alpha = .7) +
    geom_point(aes(y = Y), size = 2.6) + geom_point(aes(y = Xb), size = 2.6, shape = 21, fill = "white", colour = "#1F4E9C", stroke = 1.1) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Y, the observed SNP frequencies (black)   vs   X b, what this pie predicts (blue)",
         subtitle = "red: the residuals Y - X b", x = "SNPs in the window", y = "allele frequency") +
    theme_bw(BASE) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = BASE),
                           plot.subtitle = element_text(size = BASE - 2))
  lc <- ggplot(tibble(it = 0:NIT, L = L), aes(it, L)) + geom_line(colour = "grey50") +
    geom_point(data = tibble(it = i - 1, L = L[i]), size = 4.5, colour = "#C0392B") +
    scale_y_log10() +
    labs(title = expression(bold("loss  " * sum((Y - Xb)^2))), subtitle = sprintf("step %d", i - 1),
         x = "step", y = NULL) +
    theme_bw(BASE) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = BASE + 1),
                           plot.subtitle = element_text(size = BASE - 1))
  ((pie | lc) + plot_layout(widths = c(1, 1.3))) / fit + plot_layout(heights = c(1.1, 1)) +
    plot_annotation(title = "Find the weights b that minimise the loss",
                    subtitle = "exactly what training any model does -- here the model is eight founders mixing into a pool",
                    theme = theme(plot.title = element_text(size = BASE + 4, face = "bold"),
                                  plot.subtitle = element_text(size = BASE)))
}
dir.create("temp_aging/analysis/Talk_loss_frames", showWarnings = FALSE)
for (i in 1:(NIT + 1)) ggsave(sprintf("temp_aging/analysis/Talk_loss_frames/f%03d.png", i), frame(i),
                              width = 12, height = 8, dpi = 110, bg = "white")
imgs <- image_read(sprintf("temp_aging/analysis/Talk_loss_frames/f%03d.png", 1:(NIT + 1)))
delays <- c(rep(60, 1), rep(12, NIT - 1), 200)   # hold the first and last
gif <- image_animate(image_join(imgs), delay = delays, optimize = TRUE)
image_write(gif, "temp_aging/analysis/Talk_loss.gif")
strip <- image_append(image_scale(imgs[c(1, 16, NIT + 1)], "x600"))
image_write(strip, "temp_aging/analysis/Talk_loss_strip.png")
cat("Wrote Talk_loss.gif, Talk_loss_strip.png, Talk_loss_frames/\n")
