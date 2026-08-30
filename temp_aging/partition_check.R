# partition_check.R — can the variance partition be run on h2_rep?
#
#   Rscript temp_aging/partition_check.R
#
# RUN ON HPC3 (via reproduce.sh / run_numbers.sh). Reads the eight half-scan
# h2_falconer files, which reproduce.sh writes under scope h2 or all.
#
# The partition decomposes h2 across 8 measures (4 treatments x 2 halves) and
# subtracts the within-cell odd-even difference as a pure error term. Two
# questions have to be answered before rerunning it on h2_rep:
#
#   1. Do the halves agree well enough for the decomposition to mean anything?
#      Each half has 5 replicates, so its var/n correction is noisier. Reported
#      as the odd-even correlation, split by whether a window carries signal --
#      on a null window the halves SHOULD disagree, so a global correlation is
#      not the test.
#   2. How much of the between-cell variance is error? If the error term
#      swamps the treatment terms the partition has nothing left to report.
#
# Writes nothing; prints. run_numbers.sh captures it to numbers/partition_check.txt
# with provenance, and that file comes back through git.

suppressMessages(library(tidyverse))

SCANS <- c("AGE_SY10_F","AGE_SY20_F","AGE_SY10_M","AGE_SY20_M")
HALF  <- "process/AGE_SY_splithalf/Scans"

read_half <- function(sc, hf) {
  f <- file.path(HALF, paste0(sc, "_no89_", hf), paste0(sc, "_no89_", hf, ".h2_falconer.txt"))
  if (!file.exists(f)) return(NULL)
  read.table(f, header = TRUE, colClasses = c(sex = "character")) %>% as_tibble() %>%
    transmute(chr = CHROM, pos, sugar = str_extract(sc, "SY[12]0"),
              sex = str_sub(sc, -1), half = hf, h2 = h2_rep)
}
d <- expand_grid(sc = SCANS, hf = c("odd","even")) %>%
  pmap_dfr(function(sc, hf) read_half(sc, hf))
if (!nrow(d)) stop("no half-scan h2_falconer files; run reproduce.sh h2 first")

w <- d %>% pivot_wider(names_from = half, values_from = h2) %>%
  filter(!is.na(odd), !is.na(even)) %>%
  mutate(arm  = if_else(chr == "chrX", "chrX", "autosome"),
         mean = (odd + even)/2,
         err  = 0.5*(odd - even)^2)          # pure error, 1 df

cat(sprintf("windows x treatment: %d   (%d scans x %d windows)\n\n",
            nrow(w), n_distinct(paste(w$sugar, w$sex)), n_distinct(paste(w$chr, w$pos))))

cat("=== 1. do the halves agree? ===\n")
cat("   split by the size of the signal -- halves SHOULD disagree where there is none\n")
w %>% mutate(band = cut(mean, c(-Inf, 0.1, 0.25, 0.5, 1, Inf),
                        labels = c("<0.1","0.1-0.25","0.25-0.5","0.5-1",">1"))) %>%
  group_by(band) %>%
  summarise(n = n(), r = round(cor(odd, even), 3),
            med_err = round(median(err), 4), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n   by arm, windows carrying signal (mean h2 > 0.5):\n")
w %>% filter(mean > 0.5) %>% group_by(arm) %>%
  summarise(n = n(), r = round(cor(odd, even), 3), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n=== 2. how much of the between-cell variance is error? ===\n")
p <- w %>% group_by(chr, pos) %>% filter(n() == 4) %>%
  summarise(gm = mean(mean),
            ss_sex  = 4*(mean(mean[sex=="M"]) - mean(mean[sex=="F"]))^2/2,
            ss_diet = 4*(mean(mean[sugar=="SY20"]) - mean(mean[sugar=="SY10"]))^2/2,
            ss_err  = sum(err)/2,
            .groups = "drop")
p %>% summarise(windows = n(),
                `mean SS_sex`  = round(mean(ss_sex), 4),
                `mean SS_diet` = round(mean(ss_diet), 4),
                `mean SS_err`  = round(mean(ss_err), 4),
                `sex > err in` = paste0(round(100*mean(ss_sex > ss_err), 1), "% of windows"),
                `diet > err in`= paste0(round(100*mean(ss_diet > ss_err), 1), "% of windows")) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n   over the windows that carry signal (grand mean h2 > 0.5):\n")
p %>% filter(gm > 0.5) %>%
  summarise(windows = n(),
            `mean SS_sex` = round(mean(ss_sex), 4),
            `mean SS_diet`= round(mean(ss_diet), 4),
            `mean SS_err` = round(mean(ss_err), 4),
            `sex > err in`= paste0(round(100*mean(ss_sex > ss_err), 1), "% of windows")) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n=== 3. the fractions, error subtracted, over signal-bearing windows ===\n")
tot <- p %>% filter(gm > 0.5) %>%
  summarise(sex = mean(ss_sex) - mean(ss_err), diet = mean(ss_diet) - mean(ss_err),
            shared = mean(gm^2*4))
cat(sprintf("   shared %.3f   sex %.3f   diet %.3f   -> sex %.1f%%  diet %.1f%% of the total\n",
            tot$shared, tot$sex, tot$diet,
            100*tot$sex/(tot$shared+tot$sex+tot$diet),
            100*tot$diet/(tot$shared+tot$sex+tot$diet)))
cat("\n   (compare the current published partition: shared 85.3, sex 13.3, diet 1.7)\n")
