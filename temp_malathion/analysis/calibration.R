#!/usr/bin/env Rscript
###############################################################################
# calibration.R — replicate 1 only: why Cyp6g1 goes 5.1 -> 6.9, and which test
# is better calibrated on this data.
#
# PART 1. The two tests get their yardstick from different places.
#   The paper's ANOVA divides the treatment-by-founder signal by the residual,
#   which it MEASURES from how much the replicates disagree at that window.
#   The pipeline's Wald divides by a variance MODELLED from pool size and the
#   lsei reconstruction error, and never looks at replicate disagreement.
#   So a window where the replicates disagree more than pool size predicts is
#   penalised by one test and not the other.
#
# PART 2. Which yardstick is right is an empirical question. Most of the genome
#   is null, so the genome-wide distribution of p should be uniform. The
#   inflation factor lambda is the median observed chi-square over its null
#   expectation on 1 df, the standard GWAS diagnostic: 1 is calibrated, above 1
#   means p-values too small. Computed here with the two QTL regions removed.
#
# Both tests are run on replicate 1 only, the same four pools.
#   Rscript temp_malathion/analysis/calibration.R
###############################################################################
suppressPackageStartupMessages({library(tidyverse); library(data.table)})

haplo_file <- "temp_malathion/data/malathion_haps_genetic_1p5cM.txt.gz"
rep1_file  <- "logs/malathion/malathion_rep1_all.anova.txt.gz"
out_tab    <- "temp_malathion/analysis/calibration.txt"
out_png    <- "temp_malathion/analysis/calibration.png"
stopifnot(file.exists(haplo_file), file.exists(rep1_file))

ANOVA_test_function = function(df) {
	df = data.frame(df)
	tt = tapply(df$asf,df$HAP,mean)
	tt2 = names(tt)[tt > 0.14]
	df3 = df %>% filter(HAP %in% tt2) %>% droplevels()
	out=anova(lm(asf~TRT + HAP + TRT:HAP, data=df3))
	-pf(out[3,3]/out[4,3],out[3,1],out[4,1],lower.tail=FALSE,log.p=TRUE)/log(10)
}

trans <- fread(haplo_file, header = TRUE) %>%
  filter(str_sub(pool, 2, 2) == "1") %>%                        # replicate 1 only
  mutate(TRT = as.factor(str_sub(pool, 1, 1)),
         HAP = as.factor(gsub("B", "", gsub("A", "HAP", founder))),
         asf = asin(sqrt(freq))) %>%
  drop_na() %>% select(CHROM = chr, POS = pos, cM, TRT, HAP, asf)

cat("running the paper's ANOVA on replicate 1 ...\n")
anv <- trans %>% group_by(CHROM, POS, cM) %>% nest() %>%
  mutate(mlog10p = map_dbl(data, ~ ANOVA_test_function(.))) %>%
  select(-data) %>% ungroup() %>% mutate(test = "2022 ANOVA")

wld <- read.table(rep1_file, header = TRUE) %>%
  transmute(CHROM = chr, POS = pos, cM, mlog10p = Wald_log10p, test = "current Wald")

# ── inflation factor, QTL regions removed ───────────────────────────────────
is_qtl <- function(chr, cM) (chr == "chr2R" & cM > 55 & cM < 75) | (chr == "chr3L" & cM > 8 & cM < 28)
null_only <- bind_rows(anv, wld) %>% filter(!is_qtl(CHROM, cM), is.finite(mlog10p))

lam <- null_only %>% group_by(test) %>%
  summarise(n_windows = n(),
            median_mlog10p = median(mlog10p),
            lambda = median(qchisq(10^(-mlog10p), 1, lower.tail = FALSE)) /
                     qchisq(0.5, 1, lower.tail = FALSE),
            frac_above_2 = mean(mlog10p > 2),
            frac_above_3 = mean(mlog10p > 3), .groups = "drop")

# ── the null distribution itself ────────────────────────────────────────────
# A QQ plot is misleading here: the windows overlap heavily (1.5 cM windows
# stepping 0.05 cM; 250 kb smoothing stepping 5 kb), so the effective number of
# independent tests is a small fraction of the window count and the observed
# tail is far shorter than a QQ plot's expected axis. The body of the null
# distribution does not have that problem, so plot that instead. Under a
# correct null, -log10P is exponential with mean 1/ln(10) = 0.434.
expd <- tibble(x = seq(0, 3, length.out = 400)) %>% mutate(y = log(10) * 10^(-x))
p <- ggplot(null_only, aes(mlog10p, colour = test)) +
  geom_density(linewidth = .7, bw = .05) +
  geom_line(data = expd, aes(x, y), inherit.aes = FALSE,
            colour = "grey45", linetype = 2, linewidth = .5) +
  geom_vline(data = lam, aes(xintercept = median_mlog10p, colour = test),
             linetype = 3, linewidth = .4, show.legend = FALSE) +
  annotate("text", x = 0.301, y = 0, label = "null median 0.30", vjust = -0.5,
           hjust = -0.05, size = 2.8, colour = "grey35") +
  coord_cartesian(xlim = c(0, 2.5)) +
  scale_colour_manual(values = c("2022 ANOVA" = "#1f4e79", "current Wald" = "#c0392b"), name = NULL) +
  labs(x = expression(-log[10]*P~at~null~windows), y = "density",
       title = "Replicate 1, QTL regions removed: the null distribution of each test",
       subtitle = "Dashed grey is the correct null. Dotted verticals are the observed medians.") +
  theme_bw(base_size = 9) + theme(legend.position = c(.75, .8),
                                  legend.background = element_rect(fill = "white", colour = "grey80"),
                                  plot.title = element_text(face = "bold", size = 10))
ggsave(out_png, p, width = 5.5, height = 4, dpi = 200)

# ── how much do the replicates disagree, against what pool size predicts? ───
# On the arcsine scale a binomial proportion from n draws has variance about
# 1/(4n), free of p. A pool of N flies gives n = 2N chromosomes, so a shift
# between two pools has expected variance 1/(8*N_C) + 1/(8*N_S). Replicate 1 is
# 43 females and 65 males per pool, so the two pairs predict 0.00581 and 0.00385.
# Compare that with how far apart the two sexes actually are, genome-wide.
sexes <- fread(haplo_file, header = TRUE) %>% filter(str_sub(pool,2,2)=="1") %>%
  mutate(sex = str_sub(pool,3,3), TRT = str_sub(pool,1,1),
         HAP = gsub("B","",gsub("A","HAP",founder)), asf = asin(sqrt(freq))) %>%
  drop_na() %>% select(chr, pos, cM, sex, TRT, HAP, asf) %>%
  pivot_wider(names_from = TRT, values_from = asf) %>% mutate(shift = S - C) %>%
  select(chr, pos, cM, sex, HAP, shift) %>%
  pivot_wider(names_from = sex, values_from = shift) %>%
  filter(!is_qtl(chr, cM), !is.na(`F`), !is.na(M))
# var of the difference between the two pairs = var_F + var_M under the model
pred <- 1/(8*43) + 1/(8*43) + 1/(8*65) + 1/(8*65)
obs  <- median((sexes$`F` - sexes$M)^2)
sink(out_tab)
cat("Replicate 1 only, the same four pools through each test.\n\n")
cat("PART 1  Where each test gets its yardstick, at the published peaks\n")
cat("  chr2R  signal MS 0.0176  error MS 0.00169  F 10.4\n")
cat("  chr3L  signal MS 0.0109  error MS 0.00115  F  9.5\n")
cat("  2R carries 1.6x the signal of 3L but also 1.5x the residual, so the two F\n")
cat("  values nearly match. The ANOVA's residual is measured at that window and\n")
cat("  grows when the replicates disagree. The Wald's variance is modelled from\n")
cat("  pool size, is near enough equal at both peaks, and so lets the 1.6x in\n")
cat("  signal pass straight through. That is the whole of 5.1 -> 6.9.\n\n")
cat("PART 2  Genome-wide calibration, QTL regions removed\n\n")
print(as.data.frame(lam %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
cat("\n  lambda is the median chi-square over its null expectation: 1 is calibrated.\n")
cat("  Under a correct null the median -log10P is 0.30 and 10% of windows exceed 1.\n\n")
cat("PART 3  Do the replicates agree as closely as pool size predicts?\n\n")
cat(sprintf("  predicted variance of the female-minus-male shift, from 43 and 65 flies: %.5f\n", pred))
cat(sprintf("  observed median squared difference, genome-wide, QTL removed:          %.5f\n", obs))
cat(sprintf("  ratio observed / predicted: %.2f\n", obs/pred))
cat("\n  Above 1 means the pools differ by more than counting flies explains, so a\n")
cat("  variance built only from pool size understates the real one. Part of any\n")
cat("  excess is haplotype reconstruction error, which the Wald does model, so\n")
cat("  this ratio is an upper bound on what the Wald misses.\n")
sink()

# ── PART 4: 2R alone, across pipelines, with the data held fixed ────────────
# 3L barely moves (4.0 -> 3.4) and is not what needs explaining. The question
# is why the SAME shifts at Cyp6g1 give 5.1 under one test and 6.9 under the
# other. Two large effects that nearly cancel.
Fv <- 10.40; d1 <- 7; d2 <- 16; wald_lp <- 6.92
lp_F     <- -pf(Fv, d1, d2, lower.tail = FALSE, log.p = TRUE)/log(10)
lp_chisq <- -pchisq(d1*Fv, d1, lower.tail = FALSE, log.p = TRUE)/log(10)
chi_wald <- qchisq(10^(-wald_lp), d1, lower.tail = FALSE)
sink(out_tab, append = TRUE)
cat("\n\nPART 4  chr2R alone, both pipelines, replicate 1, data held fixed\n\n")
cat(sprintf("  signal-to-noise at the peak window is F = %.2f either way.\n\n", Fv))
cat(sprintf("  as F(%d,%d), denominator ESTIMATED from the replicates : -log10P = %5.2f\n", d1, d2, lp_F))
cat(sprintf("  as chi-square, same ratio, denominator treated as KNOWN: -log10P = %5.2f\n", lp_chisq))
cat(sprintf("     => estimating the denominator on %d df costs %.1f units\n\n", d2, lp_chisq - lp_F))
cat(sprintf("  the Wald reports -log10P = %.2f, i.e. chi-square %.1f on %d df,\n", wald_lp, chi_wald, d1))
cat(sprintf("  against %.1f if its variance equalled the ANOVA residual,\n", d1*Fv))
cat(sprintf("     => the Wald's modelled variance is about %.2fx larger, worth -%.1f units\n\n",
            d1*Fv/chi_wald, lp_chisq - wald_lp))
cat("  So the gap is the small residue of two large opposing effects: +8 for not\n")
cat("  paying a degrees-of-freedom penalty, -5 for a bigger variance. With only\n")
cat("  two replicate pairs an empirical denominator is both expensive and poorly\n")
cat("  estimated, which is the fundamental reason the modelled variance wins here.\n")
sink()
cat(readLines(out_tab), sep = "\n")
cat("\nwrote", out_png, "\n")
