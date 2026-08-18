# gather_scans.R — collapse every scan into ONE small file to bring home.
#
# RUN ON HPC3 from the repo root, once run_all.sh's jobs have finished:
#   module load R/4.2.2
#   Rscript scripts_oneoffs/AGE_2024/gather_scans.R
#
# Then fetch the single file it names:
#   scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/process/AGE_2024/AGE_2024_vs_AGE_SY.txt.gz \
#       process/AGE_2024/
#
# Nine scans go in:
#   AGE_2024                    the pilot -- 6 cages, females, lab food
#   AGE_SY{10,20}_{F,M}_R1to6   AGE_SY cut to its first 6 replicates (matched)
#   AGE_SY{10,20}_{F,M}         AGE_SY on all 12 (free, and worth having)
#
# Missing scans are reported and skipped rather than aborting, so this is usable
# while the last jobs are still running.

suppressMessages(library(tidyverse))

OUT <- "process/AGE_2024/AGE_2024_vs_AGE_SY.txt.gz"

scans <- bind_rows(
  tibble(scan = "AGE_2024", dir = "process/AGE_2024",
         diet = "lab", sex = "F", reps = "1-6"),
  expand_grid(sugar = c("SY10", "SY20"), sex = c("F", "M")) %>%
    mutate(scan = paste0("AGE_", sugar, "_", sex, "_R1to6"), dir = "process/AGE_SY",
           diet = sugar, reps = "1-6") %>% select(scan, dir, diet, sex, reps),
  expand_grid(sugar = c("SY10", "SY20"), sex = c("F", "M")) %>%
    mutate(scan = paste0("AGE_", sugar, "_", sex), dir = "process/AGE_SY",
           diet = sugar, reps = "1-12") %>% select(scan, dir, diet, sex, reps)
) %>% mutate(file = file.path(dir, "Scans", scan, paste0(scan, ".scan.txt")))

have <- scans %>% filter(file.exists(file))
gone <- scans %>% filter(!file.exists(file))
if (nrow(gone)) {
  cat("NOT FOUND (skipped -- still running?):\n")
  cat(paste0("  ", gone$file, collapse = "\n"), "\n\n")
}
if (!nrow(have)) stop("no scans found at all")

long <- have %>% pmap_dfr(function(scan, dir, diet, sex, reps, file) {
  d <- read.table(file, header = TRUE) %>% as_tibble()
  # bias columns exist only in scans made after XQTL2 #34; carry them if present
  for (v in c("Cutl_H2_bias", "Falc_H2_bias", "Cutl_clamp_frac", "cM"))
    if (!v %in% names(d)) d[[v]] <- NA_real_
  d %>% transmute(chr, pos = as.integer(pos), cM,
                  scan, diet, sex, reps,
                  Wald_log10p, Cutl_H2, Falc_H2,
                  Cutl_H2_bias, Falc_H2_bias, Cutl_clamp_frac)
})

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.table(long, gzfile(OUT), row.names = FALSE, quote = FALSE, sep = "\t")

cat("Wrote:", OUT, "\n")
cat(sprintf("  %d rows, %d windows, %d scans\n",
            nrow(long), n_distinct(paste(long$chr, long$pos)), n_distinct(long$scan)))
long %>% group_by(scan, diet, sex, reps) %>%
  summarise(windows = n(), max_wald = round(max(Wald_log10p, na.rm = TRUE), 1),
            median_h2 = round(median(Cutl_H2, na.rm = TRUE), 3),
            has_bias = !all(is.na(Cutl_H2_bias)), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
