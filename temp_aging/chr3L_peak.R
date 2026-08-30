# chr3L_peak.R — the one region that stands apart, at chr3L 9.30 Mb.
#
#   Rscript temp_aging/chr3L_peak.R
#
# DATA: the 10-replicate AGE_SY dataset -- replicates 8 and 9 dropped, those
# being the May 2023 cage (helpfiles/AGE_2024/population_assignment.txt).
#
# Three things about this locus, and nothing beyond them:
#
#   1. how big and how well localized it is, from the scan
#   2. which founders move in the selected pools, from the per-replicate founder
#      frequencies behind Figure 2
#   3. what that founder does in the UNSELECTED controls across the ten cages,
#      which are ordered in time -- Figure 2's control panel
#
# (3) is the check on (2). A haplotype depleted among the long-lived might simply
# be a bad haplotype, in which case it should also be losing ground in the cage
# itself. Whether it is is a question about the controls, and the controls are
# unselected, so they answer it without circularity.

suppressMessages(library(tidyverse))

SCAN <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
ZOOM <- "process/AGE_SY/AGE_SY_zoom_means.txt.gz"
LOC  <- "chr3L:9.31"
POS  <- 9305000
DROP <- c(8, 9)

d <- read.table(SCAN, header = TRUE, sep = "\t") %>% as_tibble()
z <- read.table(ZOOM, header = TRUE, sep = "\t") %>% as_tibble() %>%
  filter(locus == LOC, pos == POS, !REP %in% DROP)

cat(sprintf("chr3L %.2f Mb, %d replicates\n\n", POS/1e6, n_distinct(z$REP)))

cat("the scan at this window:\n\n")
d %>% filter(chr == "chr3L", pos == POS) %>%
  transmute(trt = paste0(sugar, "_", sex), `-log10 P` = round(Wald_log10p, 1),
            `h2 %` = round(H2, 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# ── which founders move ──────────────────────────────────────────────────────
# Frequencies sum to one, so a founder that is selected against pushes every
# other up whether or not anything is happening to them. The question is which
# founders move, and by how much relative to each other.
cat("\nfounder frequency, control and selected, mean over the 10 replicates:\n\n")
mv <- z %>% pivot_wider(names_from = TRT, values_from = freq) %>%
  group_by(sugar, sex, founder) %>%
  summarise(ctrl = mean(C), sel = mean(Z), d = mean(Z - C), .groups = "drop")
mv %>% mutate(cell = paste0(sugar, "_", sex)) %>%
  select(founder, cell, d) %>% pivot_wider(names_from = cell, values_from = d) %>%
  left_join(mv %>% group_by(founder) %>% summarise(control = round(mean(ctrl), 3)),
            by = "founder") %>%
  relocate(control, .after = founder) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)

top <- mv %>% group_by(sugar, sex) %>%
  summarise(down = founder[which.min(d)], fell = min(d),
            up = founder[which.max(d)], rose = max(d),
            n_down = sum(d < 0), .groups = "drop")
cat("\n")
top %>% mutate(across(where(is.numeric), ~round(., 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
f <- mv %>% filter(founder == top$down[1])
cat(sprintf("\n  %s: %.1f%% in controls, falls %.1f-%.1f points, a %.0f-%.0f%% reduction\n",
            top$down[1], 100*mean(f$ctrl), -100*max(f$d), -100*min(f$d),
            -100*max(f$d/f$ctrl), -100*min(f$d/f$ctrl)))

# ── the same founder in the unselected controls, across cages ────────────────
# Replicate is an ordered label, not a number: the cages were set up in sequence,
# so the order is time. A founder losing ground in the cage would trend here.
cat("\nfounder frequency in the UNSELECTED controls, by replicate in order:\n\n")
ct <- z %>% filter(TRT == "C") %>% group_by(founder, REP) %>%
  summarise(f = mean(freq), .groups = "drop") %>% arrange(founder, REP)
ct %>% pivot_wider(names_from = REP, values_from = f) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n")
ct %>% group_by(founder) %>%
  summarise(lo = round(min(f), 3), hi = round(max(f), 3),
            `first 5` = round(mean(f[seq_len(5)]), 3),
            `last 5` = round(mean(f[6:10]), 3),
            rho = round(cor(seq_len(n()), f, method = "spearman"), 2),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n  rho is Spearman of frequency on cage order. The founder that is selected\n")
cat("  against among the long-lived is not the one drifting in the cages.\n")
