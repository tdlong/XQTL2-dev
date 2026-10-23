# make_talk_mixing.R -- two talk slides on estimating founder frequencies from
# pool-seq, drawn as pictures rather than Y = Xb + e.
#
#   Rscript temp_aging/analysis/make_talk_mixing.R
#
# Talk_Fig_mixing.png  Y = X b as coloured blocks. Y is a column of SNP
#   frequencies (grey), X the founder genotypes with each alt allele painted in
#   its founder's colour, b the founder frequencies as a stacked bar. One SNP row
#   is pulled out to show its frequency is just the sum of the b's of the
#   founders carrying the alt allele.
# Talk_Fig_impute.png  back-imputation: raw pool-seq frequencies at 60x scatter
#   around the truth with binomial error; X b-hat, fit to all SNPs in the
#   window, lands much closer.
#
# Toy data: 8 founders, b from the chr3L 9.30 control frequencies, 16 SNPs with
# alt alleles assigned at random.

suppressMessages({library(tidyverse); library(patchwork); library(quadprog)})
source("temp_aging/analysis/palette.R")
OUT1 <- "temp_aging/analysis/Talk_Fig_mixing.png"
OUT2 <- "temp_aging/analysis/Talk_Fig_impute.png"
BASE <- 15
set.seed(7)

fnd <- names(FOUNDER_COL)
b   <- c(AB8 = .052, B1 = .145, B2 = .031, B3 = .151, B4 = .016, B5 = .041, B6 = .539, B7 = .024)
S   <- 16
X   <- matrix(rbinom(S * 8, 1, 0.4), S, 8, dimnames = list(paste0("SNP ", 1:S), fnd))
X[X %*% rep(1, 8) == 0, sample(8, 1)] <- 1              # no monomorphic rows
Y   <- as.numeric(X %*% b)
ROW <- 7                                                 # the SNP pulled out below
th  <- theme_void(BASE) + theme(plot.title = element_text(hjust = .5, face = "bold", size = BASE + 1),
                                plot.subtitle = element_text(hjust = .5, size = BASE - 1))

# ── slide 1: Y, X, and the b that fits ──────────────────────────────────────
# Y is observed with pool-seq noise, so no b reproduces it exactly: b is the
# vector that makes X b as close to Y as possible. Two candidates are shown
# side by side -- equal founder frequencies, and the least-squares b -- each
# with its X b in the same grey scale as Y, its residuals, and its sum of squares.
COV1 <- 60
Yobs <- rbinom(S, COV1, Y) / COV1
fitb <- function(y) { f <- solve.QP(Dmat = crossprod(X) + diag(1e-8, 8), dvec = crossprod(X, y),
  Amat = cbind(rep(1, 8), diag(8)), bvec = c(1, rep(0, 8)), meq = 1); setNames(f$solution, fnd) }
b_eq  <- setNames(rep(1/8, 8), fnd)
b_hat <- fitb(Yobs)
snpf  <- factor(rownames(X), levels = rev(rownames(X)))
grey  <- scale_fill_gradient(low = "grey97", high = "grey30", limits = c(0, 1), guide = "none")
col_of <- function(v, title, sub = NULL, digits = 2) ggplot(tibble(snp = snpf, v = v), aes(1, snp)) +
  geom_tile(aes(fill = v), colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf(paste0("%.", digits, "f"), v)), size = 3.4, colour = ifelse(v > .55, "white", "black")) +
  grey + coord_equal(clip = "off") + labs(title = title, subtitle = sub) + th
res_of <- function(r, title) ggplot(tibble(snp = snpf, r = r), aes(1, snp)) +
  geom_tile(aes(fill = r), colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%+.2f", r)), size = 3.2) +
  scale_fill_gradient2(low = "#1F4E9C", mid = "white", high = "#C0392B", limits = c(-.5, .5), guide = "none") +
  coord_equal(clip = "off") + labs(title = title, subtitle = sprintf("sum of squares\n%.3f", sum(r^2))) + th
bcol <- function(bv, title) ggplot(tibble(founder = factor(fnd, levels = rev(fnd)), b = bv), aes(1, founder)) +
  geom_tile(fill = "white", colour = "grey70", linewidth = .4) +
  geom_text(aes(label = sprintf("%s  %.2f", founder, b)), size = 3.4) +
  coord_equal(clip = "off") + labs(title = title) + th

Xl <- as_tibble(X, rownames = "snp") %>% pivot_longer(-snp, names_to = "founder", values_to = "alt") %>%
  mutate(snp = factor(snp, levels = rev(rownames(X))), founder = factor(founder, levels = fnd))
pX <- ggplot(Xl, aes(founder, snp)) +
  geom_tile(aes(fill = factor(alt)), colour = "grey70", linewidth = .4) +
  geom_text(aes(label = alt, colour = factor(alt)), size = 3.4) +
  scale_fill_manual(values = c(`0` = "white", `1` = "black"), guide = "none") +
  scale_colour_manual(values = c(`0` = "grey55", `1` = "white"), guide = "none") +
  scale_x_discrete(position = "top") + coord_equal() +
  labs(title = "X", subtitle = "founder genotypes, known") + th +
  theme(axis.text.x = element_text(size = BASE - 4), axis.text.y = element_text(size = BASE - 5, hjust = 1))
eq <- function(s) ggplot() + annotate("text", 0, 0, label = s, size = 11) + th
blank <- ggplot() + th

pY   <- col_of(Yobs, "Y", "observed SNP\nfrequencies, pool-seq")
left <- (pY | eq("\u2248") | pX | eq("\u00d7") | bcol(b_hat, "b, fitted")) + plot_layout(widths = c(1.3, .4, 4.2, .4, 1.6))
mid  <- (col_of(as.numeric(X %*% b_eq), "X b", "with b = 1/8\nfor every founder") | res_of(Yobs - as.numeric(X %*% b_eq), "Y \u2212 X b") |
         blank | col_of(as.numeric(X %*% b_hat), "X b", "with the\nfitted b") | res_of(Yobs - as.numeric(X %*% b_hat), "Y \u2212 X b")) +
        plot_layout(widths = c(1.3, 1.3, .5, 1.3, 1.3))
png(OUT1, width = 14, height = 7.5, units = "in", res = 250)
print(wrap_plots(wrap_elements(full = left), blank, wrap_elements(full = mid), nrow = 1, widths = c(7.9, .4, 5.7)) +
      plot_annotation(title = "Local founder frequencies from pool-seq, founders known",
                      subtitle = "b is the vector of founder frequencies, each \u2265 0 and summing to 1, that makes X b closest to Y:\nthe one that minimises the sum of squared residuals. Right: a guess of equal founders, and the fitted b.",
                      theme = theme(plot.title = element_text(size = BASE + 3, face = "bold"),
                                    plot.subtitle = element_text(size = BASE - 2),
                                    plot.margin = margin(5, 10, 5, 20))))
invisible(dev.off())
cat(sprintf("slide 1: fitted b vs true b, max abs error %.3f\n", max(abs(b_hat - b))))

# ── slide 2: back-imputation beats the raw frequencies ───────────────────────
S2 <- 60; COV <- 60
X2 <- matrix(rbinom(S2 * 8, 1, 0.4), S2, 8); X2[X2 %*% rep(1, 8) == 0, 1] <- 1
truth <- as.numeric(X2 %*% b)
raw   <- rbinom(S2, COV, truth) / COV
# constrained least squares: b >= 0, sum b = 1
fit <- solve.QP(Dmat = crossprod(X2) + diag(1e-8, 8), dvec = crossprod(X2, raw),
                Amat = cbind(rep(1, 8), diag(8)), bvec = c(1, rep(0, 8)), meq = 1)
bhat <- fit$solution; imp <- as.numeric(X2 %*% bhat)
rmse <- function(a) sqrt(mean((a - truth)^2))
d2 <- tibble(snp = 1:S2, truth, raw, imputed = imp) %>%
  pivot_longer(-c(snp, truth), names_to = "what", values_to = "freq") %>%
  mutate(what = factor(what, levels = c("raw", "imputed"),
                       labels = c(sprintf("raw pool-seq, %dx   (rmse %.3f)", COV, rmse(raw)),
                                  sprintf("imputed from the fitted founder frequencies   (rmse %.3f)", rmse(imp)))))
p2 <- ggplot(d2, aes(snp, freq)) +
  geom_line(aes(y = truth), colour = "grey55", linewidth = .7) +
  geom_point(aes(colour = what, shape = what), size = 2.4) +
  scale_colour_manual(values = c("#D55E00", "#1F4E9C"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "SNPs along the window", y = "alt allele frequency",
       title = "Back-imputed SNP frequencies are closer to the truth than the raw data",
       subtitle = sprintf("grey line: true frequency. %d SNPs estimate 8 founder frequencies, so each imputed value borrows strength from all %d.", S2, S2)) +
  theme_bw(BASE) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                         legend.text = element_text(size = BASE - 1))
png(OUT2, width = 13, height = 6.5, units = "in", res = 250)
print(p2)
invisible(dev.off())
cat("Wrote:", OUT1, "and", OUT2, "\n")
cat(sprintf("slide 2: rmse raw %.3f, imputed %.3f; b-hat vs b max abs error %.3f\n", rmse(raw), rmse(imp), max(abs(bhat - b))))
