# gather_splithalf_H2.R — the eight split-half scans into one long data frame.
#
# Output: process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz
#   chr  pos  sex  sugar  half  Cutl_H2  Falc_H2  Wald_log10p
#
# One row per window per cell: 5 chromosomes x (2 sexes x 2 sugars x 2 halves).
# Falc_H2 and Wald_log10p ride along because they cost nothing to carry and
# save a second trip if you want to check the halves against the main scan.
#
# Small enough to pull down and explore locally -- that is the point of it.
#
# Run from the XQTL2-dev repo ROOT (the submit script chains this automatically):
#   module load R/4.2.2
#   Rscript scripts_oneoffs/AGE_SY/splithalf/gather_splithalf_H2.R

suppressMessages(library(tidyverse))

BASE <- "process/AGE_SY_splithalf/Scans"
OUT  <- "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"

# Cells spelled out rather than parsed back off the folder name.
cells <- expand_grid(
  sugar = c("SY10", "SY20"),
  sex   = c("F", "M"),
  half  = c("odd", "even")
) %>%
  mutate(scan = paste0("AGE_", sugar, "_", sex, "_", half),
         file = file.path(BASE, scan, paste0(scan, ".scan.txt")))

missing <- cells$file[!file.exists(cells$file)]
if (length(missing))
  stop("missing scan output:\n  ", paste(missing, collapse = "\n  "))

long <- pmap_dfr(cells, function(sugar, sex, half, scan, file) {
  d <- read.table(file, header = TRUE) %>% as_tibble()
  # XQTL2 #34 added the squaring-bias columns. Carry them if present -- H2 has
  # a floor and the partition's main term rests on it -- but do not require
  # them, so scans made before #34 still gather.
  for (v in c("Cutl_H2_bias", "Falc_H2_bias", "Cutl_clamp_frac"))
    if (!v %in% names(d)) d[[v]] <- NA_real_
  d %>% transmute(chr, pos = as.integer(pos), sex, sugar, half,
                  Cutl_H2, Falc_H2, Wald_log10p,
                  Cutl_H2_bias, Falc_H2_bias, Cutl_clamp_frac)
})

if (all(is.na(long$Cutl_H2_bias)))
  cat("NOTE: no Cutl_H2_bias column -- these scans predate XQTL2 #34.\n",
      "      Re-run them to get the H2 floor.\n", sep = "")

# Every cell must cover the same windows, or the contrasts below are comparing
# different places. Fail here rather than silently recycling downstream.
per_cell <- long %>% count(sex, sugar, half, name = "n_windows")
if (n_distinct(per_cell$n_windows) != 1) {
  print(per_cell)
  stop("cells cover different numbers of windows")
}
grid_ok <- long %>%
  group_by(sex, sugar, half) %>%
  summarise(key = paste0(chr, ":", pos, collapse = ","), .groups = "drop") %>%
  pull(key) %>% n_distinct()
if (grid_ok != 1) stop("cells cover different windows (same count, different positions)")

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.table(long, gzfile(OUT), row.names = FALSE, quote = FALSE, sep = "\t")

cat("Wrote:", OUT, "\n")
cat(sprintf("  %d rows = %d windows x %d cells\n",
            nrow(long), per_cell$n_windows[1], nrow(per_cell)))
print(per_cell)
cat("\nCutl_H2 by cell (bias = the reported floor, XQTL2 #34):\n")
long %>%
  group_by(sex, sugar, half) %>%
  summarise(median = median(Cutl_H2, na.rm = TRUE),
            median_bias = median(Cutl_H2_bias, na.rm = TRUE),
            median_corrected = median(Cutl_H2 - Cutl_H2_bias, na.rm = TRUE),
            clamp_frac = median(Cutl_clamp_frac, na.rm = TRUE),
            frac_negative = mean(Cutl_H2 < 0, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
