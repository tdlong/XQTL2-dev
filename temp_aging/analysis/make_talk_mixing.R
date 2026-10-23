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

# ── slide 1: the three blocks ────────────────────────────────────────────────
Xl <- as_tibble(X, rownames = "snp") %>% pivot_longer(-snp, names_to = "founder", values_to = "alt") %>%
  mutate(snp = factor(snp, levels = rev(rownames(X))), founder = factor(founder, levels = fnd))
pX <- ggplot(Xl, aes(founder, snp)) +
  geom_tile(aes(fill = factor(alt)), colour = "grey70", linewidth = .4) +
  geom_text(aes(label = alt, colour = factor(alt)), size = 3.6) +
  geom_tile(data = Xl %>% filter(snp == paste0("SNP ", ROW)), fill = NA, colour = "red", linewidth = 1.3) +
  scale_fill_manual(values = c(`0` = "white", `1` = "black"), guide = "none") +
  scale_colour_manual(values = c(`0` = "grey55", `1` = "white"), guide = "none") +
  scale_x_discrete(position = "top") + coord_equal() +
  labs(title = "X", subtitle = "founder genotypes at each SNP, known") + th +
  theme(axis.text.x = element_text(size = BASE - 3), axis.text.y = element_text(size = BASE - 4, hjust = 1))

Yl <- tibble(snp = factor(rownames(X), levels = rev(rownames(X))), y = Y)
pY <- ggplot(Yl, aes(1, snp)) +
  geom_tile(fill = "white", colour = "grey70", linewidth = .4) +
  geom_text(aes(label = sprintf("%.2f", y)), size = 3.8) +
  geom_tile(data = Yl %>% filter(snp == paste0("SNP ", ROW)), fill = NA, colour = "red", linewidth = 1.3) +
  coord_equal(clip = "off") + labs(title = "Y", subtitle = "SNP frequencies\nin the pool, observed") + th

bl <- tibble(founder = factor(fnd, levels = rev(fnd)), b = b[fnd])
pb <- ggplot(bl, aes(1, founder)) +
  geom_tile(fill = "white", colour = "grey70", linewidth = .4) +
  geom_text(aes(label = sprintf("%s  %.2f", founder, b)), size = 3.8) +
  coord_equal(ratio = 1, clip = "off") +
  labs(title = "b", subtitle = "founder frequencies,\nestimated; sum to 1") + th

eq <- function(s) ggplot() + annotate("text", 0, 0, label = s, size = 12) + th

# the boxed row as a dot product: X row, b underneath, the products, the sum
xr <- X[ROW, ]; pr <- xr * b
strip <- bind_rows(
  tibble(row = "row of X", founder = fnd, lab = as.character(xr), on = xr == 1),
  tibble(row = "b",        founder = fnd, lab = sprintf("%.2f", b), on = xr == 1),
  tibble(row = "product",  founder = fnd, lab = ifelse(xr == 1, sprintf("%.2f", pr), "0"), on = xr == 1)) %>%
  mutate(row = factor(row, levels = rev(c("row of X", "b", "product"))), founder = factor(founder, levels = fnd))
pS <- ggplot(strip, aes(founder, row)) +
  geom_tile(fill = "white", colour = "grey70", linewidth = .4) +
  geom_text(aes(label = lab, colour = on, fontface = ifelse(on, "bold", "plain")), size = 4.2) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"), guide = "none") +
  scale_x_discrete(position = "top") +
  annotate("text", x = 8.9, y = 1, label = sprintf("sum = %.2f = Y for SNP %d", Y[ROW], ROW), hjust = 0, size = 4.6, fontface = "bold") +
  annotate("text", x = 8.9, y = 3, label = sprintf("SNP %d", ROW), hjust = 0, size = 4.2, colour = "red") +
  coord_cartesian(xlim = c(.5, 11.5), clip = "off") +
  labs(subtitle = "each SNP frequency is the sum of the b of the founders carrying its alt allele") +
  th + theme(axis.text.x = element_text(size = BASE - 4), axis.text.y = element_text(size = BASE - 3, hjust = 1),
             plot.subtitle = element_text(hjust = 0))

top <- (pY | eq("=") | pX | eq("×") | pb) + plot_layout(widths = c(1.4, .5, 5, .5, 1.8))
png(OUT1, width = 11, height = 8.5, units = "in", res = 250)
print((top / pS) + plot_layout(heights = c(6, 1.6)) +
      plot_annotation(title = "Local founder frequencies from pool-seq, founders known",
                      subtitle = "find the b that best explains Y, with b \u2265 0 and summing to 1",
                      theme = theme(plot.title = element_text(size = BASE + 3, face = "bold"),
                                    plot.subtitle = element_text(size = BASE),
                                    plot.margin = margin(5, 10, 5, 25))))
invisible(dev.off())

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
