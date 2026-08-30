# diagnose_clamp.R — does the penetrance clamp bind on founders that carry h2?
#
#   Rscript temp_aging/diagnose_clamp.R [means_file]
#
# Cutler h2 = 200 * sum_f C_f * Affect_f^2, with Pen_f = Z_f*P/C_f clamped to
# [P/2, 2P]. A founder at the lsei floor clamps constantly and contributes
# ~200*3e-4*0.36^2 = 0.008, i.e. nothing. So a high clamp FRACTION is not
# evidence of anything -- Cutler_clamp_frac counts founders unweighted.
#
# The question is whether the clamp binds on founders with real frequency, where
# one clamped founder at C=0.10 carries ~2.6 of h2. This weights the clamp by
# what each founder actually contributes.
#
# Default input is the local zoom_means (7 peaks, 1.2 Mb each, includes
# chrX:10.09). Pass a meansBySample.txt to run it genome-wide on HPC3.

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
IN   <- if (length(args) >= 1) args[1] else "process/AGE_SY/AGE_SY_zoom_means.txt.gz"
stopifnot(file.exists(IN))

# "F" must not become FALSE.
m <- read.table(IN, header = TRUE, sep = "\t", colClasses = c(sex = "character")) %>%
  as_tibble() %>% mutate(sex = ifelse(sex %in% c("FALSE", "F"), "F", "M"))

# Selected proportion per replicate, from the designs.
des <- list.files("helpfiles/AGE_SY/nov_only", "^AGE_SY[12]0_[FM]\\.no89\\.txt$",
                  full.names = TRUE) %>%
  map_dfr(~ read.table(.x, header = TRUE, colClasses = c(Sex = "character")) %>%
            as_tibble()) %>%
  filter(TRT == "Z") %>%
  transmute(sugar = longTRT, sex = ifelse(Sex %in% c("FALSE","F"), "F", "M"),
            REP = as.integer(REP), P = Proportion) %>%
  mutate(sugar = str_replace(sugar, "Age", "")) %>% distinct()

d <- m %>%
  select(chr, pos, TRT, REP, founder, freq, sugar, sex) %>%
  pivot_wider(names_from = TRT, values_from = freq) %>%
  filter(!is.na(C), !is.na(Z)) %>%
  inner_join(des, by = c("sugar", "sex", "REP")) %>%
  mutate(arm     = ifelse(chr == "chrX", "chrX", "autosome"),
         Pen_raw = Z * P / C,
         clamped = Pen_raw < P/2 | Pen_raw > 2*P,
         Pen_hat = pmin(pmax(Pen_raw, P/2), 2*P),
         Affect  = qnorm(1 - P) - qnorm(1 - Pen_hat),
         h2_con  = 200 * C * Affect^2)          # this founder's share of h2

cat("rows:", nrow(d), " windows:", n_distinct(paste(d$chr, d$pos)), "\n\n")

cat("=== clamp rate, unweighted vs weighted by h2 carried ===\n")
d %>% group_by(arm, sex) %>%
  summarise(founders      = n(),
            clamp_pct     = round(100*mean(clamped), 1),
            h2_total      = round(sum(h2_con), 0),
            h2_pct_clamped= round(100*sum(h2_con[clamped])/sum(h2_con), 1),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n=== split by whether the founder is at the floor (C < 0.01) ===\n")
d %>% mutate(cls = ifelse(C < 0.01, "at floor", "real freq")) %>%
  group_by(arm, sex, cls) %>%
  summarise(n = n(), clamp_pct = round(100*mean(clamped), 1),
            mean_h2_each = round(mean(h2_con), 3),
            h2_share = round(sum(h2_con), 0), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n=== among founders with REAL frequency (C >= 0.05), who clamps ===\n")
d %>% filter(C >= 0.05) %>% group_by(arm, sex) %>%
  summarise(n = n(), clamp_pct = round(100*mean(clamped), 1),
            med_C = round(median(C), 3),
            med_ratio = round(median(Z/C), 2),
            iqr_ratio = round(IQR(Z/C), 2), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
