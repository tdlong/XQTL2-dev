# make_4scan_df.R — collapse the four 12-replicate scans into ONE small file.
#
# RUN THIS ON HPC3, from /dfs7/adl/tdlong/fly_pool/XQTL2-dev:
#   module load R/4.2.2
#   Rscript temp_aging/make_4scan_df.R
#
# Then fetch the single file it writes:
#   scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/process/AGE_SY/AGE_SY_4scan.txt.gz \
#       process/AGE_SY/
#
# The point is to avoid scp -r on the Scans directories, which drags down the
# smoothed .rds files and the per-chromosome intermediates for the sake of two
# columns.

suppressMessages(library(tidyverse))

SCANDIR <- "process/AGE_SY/Scans"
OUT     <- "process/AGE_SY/AGE_SY_4scan.txt.gz"

cells <- expand_grid(sugar = c("SY10", "SY20"), sex = c("F", "M")) %>%
  mutate(scan = paste0("AGE_", sugar, "_", sex),
         file = file.path(SCANDIR, scan, paste0(scan, ".scan.txt")))

missing <- cells$file[!file.exists(cells$file)]
if (length(missing)) stop("missing:\n  ", paste(missing, collapse = "\n  "))

long <- pmap_dfr(cells, function(sugar, sex, scan, file) {
  d <- read.table(file, header = TRUE) %>% as_tibble()
  # bias columns exist only in scans made after XQTL2 #34; carry them if present
  for (v in c("Cutl_H2_bias", "Falc_H2_bias", "Cutl_clamp_frac"))
    if (!v %in% names(d)) d[[v]] <- NA_real_
  # cM comes from add_genetic() in the pipeline. Needed because clumping in
  # physical distance is wrong here: h2 is smeared over ~9 Mb in the
  # low-recombination blocks flanking the centromeres and over ~1 Mb in
  # euchromatin, so a fixed Mb window counts one block as ten loci.
  if (!"cM" %in% names(d)) d$cM <- NA_real_
  d %>% transmute(chr, pos = as.integer(pos), cM, sugar, sex,
                  Wald_log10p, Cutl_H2, Falc_H2,
                  Cutl_H2_bias, Falc_H2_bias, Cutl_clamp_frac)
})

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.table(long, gzfile(OUT), row.names = FALSE, quote = FALSE, sep = "\t")

cat("Wrote:", OUT, "\n")
cat(sprintf("  %d rows = %d windows x %d treatments\n",
            nrow(long), n_distinct(paste(long$chr, long$pos)), nrow(cells)))
long %>% group_by(sugar, sex) %>%
  summarise(max_wald = round(max(Wald_log10p, na.rm = TRUE), 1),
            median_h2 = round(median(Cutl_H2, na.rm = TRUE), 3),
            has_bias = !all(is.na(Cutl_H2_bias)), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
