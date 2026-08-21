# control_drift.R — is there a CONSISTENT pattern across the seven regions in
# the control frequencies, or just drift?
#
#   Rscript temp_aging/control_drift.R
#
# DATA: the 10-replicate AGE_SY dataset -- replicates 8 and 9 dropped, those
# being the May 2023 cage (helpfiles/AGE_2024/population_assignment.txt).
#
# THE QUESTION. Panel 8 of Figure 2 shows, for each of seven peaks, the haplotype
# most enriched in the long-lived pool and the one most depleted, plotted as
# their frequency in the UNSELECTED controls across replicates. The worry it
# addresses is that a haplotype might simply be leaving the cages, which would
# make it look depleted in any selected pool.
#
# WHAT IS NOT THE QUESTION. Whether a single series has a slope different from
# zero. Under drift a haplotype goes up or down; there is no expectation that any
# one of them stays flat, so a per-series test answers nothing. An earlier
# version of this file did exactly that -- fourteen regressions, then each
# compared against a simulated random walk -- and it was measuring drift against
# drift.
#
# WHAT IS. Whether the seven loci agree. If the hits were founders on their way
# out of the cages, the susceptible set would slope down ACROSS loci and the
# protective set up, consistently. Seven paired observations, one per locus.

suppressMessages(library(tidyverse))

MEANS <- "process/AGE_SY/AGE_SY_zoom_means.txt.gz"
SCAN  <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
DROP_REP <- c(8, 9)          # the zoom means carry every replicate
RARE  <- 0.025               # founders below this in the controls are noise

m <- read.table(MEANS, header = TRUE, sep = "\t") %>% as_tibble() %>%
  filter(!REP %in% DROP_REP)
s <- read.table(SCAN, header = TRUE, sep = "\t") %>% as_tibble()

# the treatment with the strongest signal at each peak, as Figure 2 draws
best <- m %>% distinct(locus, chr, peak_pos) %>%
  pmap_dfr(function(locus, chr, peak_pos)
    s %>% filter(chr == !!chr, pos == !!peak_pos) %>%
      slice_max(Wald_log10p, n = 1) %>% transmute(locus, sugar, sex))

w <- m %>% inner_join(best, by = c("locus", "sugar", "sex")) %>%
  filter(pos == peak_pos) %>%
  pivot_wider(names_from = TRT, values_from = freq, names_prefix = "f")

# the two movers at each peak, rare founders excluded
mv <- w %>% group_by(locus, founder) %>%
  summarise(D = mean(fZ - fC), mC = mean(fC), .groups = "drop") %>%
  filter(mC >= RARE) %>% group_by(locus) %>%
  { bind_rows(slice_max(., D, n = 1) %>% mutate(dir = "protective"),
              slice_min(., D, n = 1) %>% mutate(dir = "susceptible")) } %>%
  ungroup() %>% select(locus, founder, dir, D)

# slope of each mover's CONTROL frequency on replicate order
slope <- w %>% inner_join(mv, by = c("locus", "founder")) %>%
  mutate(ri = match(REP, sort(unique(REP)))) %>%
  group_by(locus, founder, dir) %>%
  summarise(slope = coef(lm(fC ~ ri))[2], .groups = "drop")

cat("per locus: the two movers and the slope of their control frequency\n")
cat("(slope is per replicate, in the order the replicates ran)\n\n")
slope %>% left_join(mv %>% select(locus, founder, D), by = c("locus","founder")) %>%
  transmute(locus, dir, founder, `change under selection` = round(D, 3),
            `control slope` = round(slope, 4)) %>%
  arrange(dir, locus) %>% as.data.frame() %>% print(row.names = FALSE)

# ── the only test: do the seven loci agree? ──────────────────────────────────
p <- slope %>% select(locus, dir, slope) %>%
  pivot_wider(names_from = dir, values_from = slope) %>%
  mutate(diff = protective - susceptible)

cat("\nIS THERE A CONSISTENT PATTERN ACROSS THE SEVEN LOCI?\n\n")
cat(sprintf("  protective slopes: %d of %d positive   (mean %+.4f)\n",
            sum(p$protective > 0), nrow(p), mean(p$protective)))
cat(sprintf("  susceptible slopes: %d of %d negative  (mean %+.4f)\n",
            sum(p$susceptible < 0), nrow(p), mean(p$susceptible)))
cat(sprintf("  protective above susceptible at %d of %d loci\n",
            sum(p$diff > 0), nrow(p)))

tt <- t.test(p$protective, p$susceptible, paired = TRUE)
bt <- binom.test(sum(p$diff > 0), nrow(p), 0.5)
cat(sprintf("\n  paired t across loci:  difference %+.4f  [%+.4f, %+.4f]  P = %.3f\n",
            tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value))
cat(sprintf("  sign test across loci: P = %.3f\n", bt$p.value))
cat("\n  With seven loci neither test has much power, so read the counts above\n")
cat("  rather than the P values: a real purging effect would show as 7 of 7,\n")
cat("  not as a marginal average.\n")
