#!/usr/bin/env Rscript
###############################################################################
# denominator_noise.R — how much of the 2022 scan's window-to-window movement
# is the denominator being re-estimated, badly, at every window?
#
# The 2022 scan estimates its error mean square separately at each window from
# 16 residual df, with no smoothing. Adjacent windows are 1.5 cM wide stepping
# 0.05 cM, so they share about 97% of their data: the TRUE error variance is
# essentially identical between neighbours. Any difference between neighbouring
# denominators is therefore estimation noise, and it goes straight into
# -log10P, raising some windows and lowering others.
#
# Three things are measured, replicate 1, published 2022 frequencies:
#   1. how variable the denominator is between neighbouring windows, against
#      what a chi-square on 16 df predicts if it were pure noise;
#   2. how much of the scan's roughness disappears if the denominator is
#      replaced by its local median and the numerator left alone;
#   3. what that does to the chr2R peak.
#
#   Rscript temp_malathion/analysis/denominator_noise.R
###############################################################################
suppressPackageStartupMessages({library(tidyverse); library(data.table)})

haplo_file <- "temp_malathion/data/malathion_haps_genetic_1p5cM.txt.gz"
out_tab <- "temp_malathion/analysis/denominator_noise.txt"
out_png <- "temp_malathion/analysis/denominator_noise.png"
stopifnot(file.exists(haplo_file))

# per window: the two mean squares the paper's F is built from
components <- function(df) {
  df <- data.frame(df)
  tt <- tapply(df$asf, df$HAP, mean); keep <- names(tt)[tt > 0.14]
  d3 <- df %>% filter(HAP %in% keep) %>% droplevels()
  a <- anova(lm(asf ~ TRT + HAP + TRT:HAP, data = d3))
  tibble(MS_signal = a[3,3], df1 = a[3,1], MS_error = a[4,3], df2 = a[4,1])
}

trans <- fread(haplo_file, header = TRUE) %>%
  filter(str_sub(pool, 2, 2) == "1") %>%
  mutate(TRT = as.factor(str_sub(pool,1,1)),
         HAP = as.factor(gsub("B","",gsub("A","HAP",founder))),
         asf = asin(sqrt(freq))) %>%
  drop_na() %>% select(CHROM = chr, POS = pos, cM, TRT, HAP, asf)

cat("decomposing the 2022 scan into its two mean squares ...\n")
mq <- trans %>% group_by(CHROM, POS, cM) %>% nest() %>%
  mutate(c = map(data, components)) %>% select(-data) %>% unnest(c) %>% ungroup() %>%
  arrange(CHROM, cM) %>%
  mutate(F_raw = MS_signal/MS_error,
         lp_raw = -pf(F_raw, df1, df2, lower.tail = FALSE, log.p = TRUE)/log(10))

# denominator replaced by a local median over +/-20 windows (+/-1 cM), numerator untouched
k <- 200          # +/-10 cM: wide enough that the averaged windows are independent
mq <- mq %>% group_by(CHROM) %>%
  mutate(MS_error_smooth = zoo::rollapply(MS_error, 2*k+1, median, partial = TRUE, align = "center"),
         F_sm  = MS_signal/MS_error_smooth,
         lp_sm = -pf(F_sm, df1, df2, lower.tail = FALSE, log.p = TRUE)/log(10)) %>% ungroup()

# 1. The denominator's noise is SHARED between neighbours, because neighbouring
# windows share 97% of their data. The scale that matters is the one at which
# windows stop sharing data: 1.5 cM, i.e. 30 steps. Measure how far the
# denominator wanders at a range of lags, against pure 16 df noise.
pred_sd <- sqrt(2) * sqrt(2/16)   # sd of log(MS1/MS2) for two independent 16 df MS
lags <- c(1, 10, 30, 60, 200)
nbr <- map_dfr(lags, function(L)
  mq %>% group_by(CHROM) %>%
    summarise(lag_steps = L, lag_cM = L * 0.05,
              sd_log_diff = sd(log(MS_error) - lag(log(MS_error), L), na.rm = TRUE),
              .groups = "drop")) %>%
  group_by(lag_steps, lag_cM) %>%
  summarise(sd_log_diff = mean(sd_log_diff), frac_of_pure_noise = sd_log_diff/pred_sd,
            .groups = "drop")

# 2. roughness of the scan, before and after
rough <- mq %>% group_by(CHROM) %>%
  summarise(raw = median(abs(diff(lp_raw))), smoothed = median(abs(diff(lp_sm))), .groups = "drop")

# 3. the chr2R peak.  A denominator borrowed from +/-10 cM is built from about
# 400/30 = 13 independent 16 df estimates, so it carries roughly 213 df, not 16.
# Judging it against F(7,16) would be the wrong reference, so report both.
df2_borrowed <- round((2*k+1)/30) * 16
mq <- mq %>% mutate(lp_sm_rightdf = -pf(F_sm, df1, df2_borrowed, lower.tail = FALSE, log.p = TRUE)/log(10))
pk <- mq %>% filter(CHROM == "chr2R", cM > 55, cM < 75) %>%
  summarise(raw_max = max(lp_raw), raw_at = cM[which.max(lp_raw)],
            sm_max = max(lp_sm),  sm_at  = cM[which.max(lp_sm)],
            sm_rightdf_max = max(lp_sm_rightdf), sm_rightdf_at = cM[which.max(lp_sm_rightdf)])

p <- mq %>% filter(CHROM == "chr2R", cM > 55, cM < 75) %>%
  select(cM, `denominator re-estimated every window (2022 scan)` = lp_raw,
         `denominator replaced by its local median` = lp_sm) %>%
  pivot_longer(-cM) %>%
  ggplot(aes(cM, value, colour = name)) + geom_line(linewidth = .5) +
  geom_vline(xintercept = 65.607, linetype = 2, colour = "grey70") +
  scale_colour_manual(values = c("#c0392b", "#1f4e79"), name = NULL) +
  labs(x = "cM on chr2R", y = expression(-log[10]*P),
       title = "The 2022 scan at Cyp6g1, numerator unchanged",
       subtitle = "Only the denominator differs: re-estimated at every window, or taken as a local median.") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom",
                                  plot.title = element_text(face = "bold", size = 10))
ggsave(out_png, p, width = 6.5, height = 4, dpi = 200)

sink(out_tab)
cat("The 2022 scan, replicate 1. Numerator and denominator separated.\n\n")
cat("1. How far does the denominator wander, and is it noise?\n")
cat("   Windows are 1.5 cM wide stepping 0.05 cM, so neighbours share ~97% of\n")
cat("   their data and their denominators are nearly the same estimate. Windows\n")
cat("   30 steps apart share nothing. If the wandering is pure estimation noise\n")
cat(sprintf("   on 16 df, sd of the log ratio reaches %.3f once windows are independent.\n\n", pred_sd))
print(as.data.frame(nbr %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\n2. Median jump in -log10P between neighbouring windows, denominator as is\n")
cat("   against the denominator replaced by a median over +/-10 cM:\n\n")
print(as.data.frame(rough %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\n3. The chr2R peak:\n\n")
print(as.data.frame(pk %>% mutate(across(everything(), ~ round(.x, 2)))), row.names = FALSE)
cat(sprintf("\n   raw_max      denominator from this window alone, 16 df, judged as F(7,16)\n"))
cat(sprintf("   sm_max       denominator borrowed from +/-10 cM, still judged as F(7,16)\n"))
cat(sprintf("   sm_rightdf   the same, judged as F(7,%d), which is the correct reference\n", df2_borrowed))
cat("\nThe numerator is identical in all three. Everything that changes is where\n")
cat("the denominator comes from, and how many df it is credited with.\n")
sink()
cat(readLines(out_tab), sep = "\n")
