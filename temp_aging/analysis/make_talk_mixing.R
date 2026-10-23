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
  mutate(snp = factor(snp, levels = rev(rownames(X))), founder = factor(founder, levels = fnd),
         fill = ifelse(alt == 1, as.character(founder), "ref"))
pX <- ggplot(Xl, aes(founder, snp)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = 1) +
  geom_tile(data = Xl %>% filter(snp == paste0("SNP ", ROW)), fill = NA, colour = "black", linewidth = 1.2) +
  scale_fill_manual(values = c(FOUNDER_COL, ref = "grey92"), guide = "none") +
  scale_x_discrete(position = "top") + coord_equal() +
  labs(title = "X", subtitle = "founder genotypes, known") + th +
  theme(axis.text.x = element_text(size = BASE - 3), axis.text.y = element_text(size = BASE - 4, hjust = 1))

Yl <- tibble(snp = factor(rownames(X), levels = rev(rownames(X))), y = Y)
pY <- ggplot(Yl, aes(1, snp)) +
  geom_tile(aes(fill = y), colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", y)), size = 3.6, colour = ifelse(Yl$y > .5, "white", "black")) +
  geom_tile(data = Yl %>% filter(snp == paste0("SNP ", ROW)), fill = NA, colour = "black", linewidth = 1.2) +
  scale_fill_gradient(low = "grey95", high = "grey25", limits = c(0, 1), guide = "none") +
  coord_equal(clip = "off") + labs(title = "Y", subtitle = "SNP frequencies\nin the pool") + th

bl <- tibble(founder = factor(fnd, levels = rev(fnd)), b = b[fnd])
pb <- ggplot(bl, aes(1, b, fill = founder)) +
  geom_col(width = .9, colour = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(b > .04, sprintf("%s %.2f", founder, b), "")),
            position = position_stack(vjust = .5), size = 3.4, colour = ifelse(rev(fnd) == "B6", "white", "black")) +
  scale_fill_manual(values = FOUNDER_COL, guide = "none") +
  coord_cartesian(xlim = c(.3, 1.7)) +
  labs(title = "b", subtitle = "founder frequencies,\nestimated; sum to 1") + th

eq <- function(s) ggplot() + annotate("text", 0, 0, label = s, size = 12) + th
carriers <- fnd[X[ROW, ] == 1]
expl <- sprintf("SNP %d is the alt allele in %s,\nso its frequency is %s = %.2f",
                ROW, paste(carriers, collapse = ", "),
                paste(sprintf("%.2f", b[carriers]), collapse = " + "), Y[ROW])
pE <- ggplot() + annotate("text", 0, 0, label = expl, size = 5, lineheight = .95) + th

top <- (pY | eq("=") | pX | eq("×") | pb) + plot_layout(widths = c(1.4, .5, 5, .5, 1.6))
png(OUT1, width = 11, height = 7.5, units = "in", res = 250)
print((top / pE) + plot_layout(heights = c(6, 1)) +
      plot_annotation(title = "Local founder frequencies from pool-seq, founders known",
                      theme = theme(plot.title = element_text(size = BASE + 3, face = "bold"),
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
