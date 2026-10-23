# make_talk_mixing.R -- candidate talk slides for "founder frequencies from
# pool-seq, founders known". Three pictures of the same idea; keep the one that
# lands in thirty seconds.
#
#   Rscript temp_aging/analysis/make_talk_mixing.R
#
#   Talk_mix_A.png  dots, line, pie: observed SNP frequencies as dots; a wrong
#                   pie draws a line that misses them, the fitted pie draws one
#                   through them.
#   Two other versions (the pool as a stack of founder chromosomes; one SNP at a
#   time with pie slices lit) were drawn and dropped -- see git history before
#   2026-10-23 if wanted.
#
# Toy data; b is the chr3L 9.30 control frequencies.

suppressMessages({library(tidyverse); library(patchwork); library(quadprog)})
source("temp_aging/analysis/palette.R")
BASE <- 16; set.seed(11)
fnd <- names(FOUNDER_COL)
b   <- c(AB8 = .052, B1 = .145, B2 = .031, B3 = .151, B4 = .016, B5 = .041, B6 = .539, B7 = .024)
S <- 40; COV <- 80
X <- matrix(rbinom(S * 8, 1, 0.4), S, 8, dimnames = list(NULL, fnd)); X[rowSums(X) == 0, 1] <- 1
truth <- as.numeric(X %*% b); Yobs <- rbinom(S, COV, truth) / COV
fitb <- function(y) { f <- solve.QP(crossprod(X) + diag(1e-8, 8), crossprod(X, y),
  cbind(rep(1, 8), diag(8)), c(1, rep(0, 8)), meq = 1); setNames(pmax(f$solution, 0), fnd) }
b_hat <- fitb(Yobs); b_eq <- setNames(rep(1/8, 8), fnd)
pie <- function(bv, title) ggplot(tibble(f = factor(fnd, levels = fnd), b = bv), aes(1, b, fill = f)) +
  geom_col(width = 1, colour = "white", linewidth = .6) + coord_polar(theta = "y") +
  scale_fill_manual(values = FOUNDER_COL, name = NULL) + labs(title = title) + theme_void(BASE) +
  theme(plot.title = element_text(hjust = .5, size = BASE), legend.position = "none")
leg <- ggplot(tibble(f = factor(fnd, levels = fnd), x = 1:8), aes(x, 1, colour = f)) + geom_point(size = 5, alpha = 0) +
  scale_colour_manual(values = FOUNDER_COL, name = "founders") + theme_void(BASE) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 5))) + theme(legend.position = "bottom")

# ── A: dots, line, pie ───────────────────────────────────────────────────────
dl <- function(bv, title, col) ggplot(tibble(snp = 1:S, obs = Yobs, fit = as.numeric(X %*% bv)), aes(snp)) +
  geom_line(aes(y = fit), colour = col, linewidth = 1.1) +
  geom_point(aes(y = obs), size = 2.6) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  labs(title = title, x = "SNPs along a 500 kb window", y = "allele frequency\nin the pool") +
  theme_bw(BASE) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = BASE))
pA <- ((dl(b_eq, "guess: every founder 1/8  —  the line misses the dots", "#C0392B") | pie(b_eq, "")) +
         plot_layout(widths = c(4, 1))) /
      ((dl(b_hat, "fitted founder frequencies  —  the line goes through them", "#1F4E9C") | pie(b_hat, "")) +
         plot_layout(widths = c(4, 1))) / leg + plot_layout(heights = c(1, 1, .15)) +
  plot_annotation(title = "Dots: SNP frequencies from sequencing the pool.  Pie: what mix of the eight founders the pool is.",
                  subtitle = "Each founder's alleles are known, so any pie predicts every dot. The founder frequencies are the pie whose line fits the dots.",
                  theme = theme(plot.title = element_text(size = BASE + 2, face = "bold"), plot.subtitle = element_text(size = BASE)))
ggsave("temp_aging/analysis/Talk_mix_A.png", pA, width = 13, height = 8, dpi = 220, bg = "white")

cat("Wrote Talk_mix_A.png\n")
