# gather_snp.R — the four SNP scans into one file.
#
# RUN ON HPC3 from the repo root, once run_snp_scans.sh has finished:
#   module load R/4.2.2
#   Rscript scripts_oneoffs/AGE_SY/nov_only/gather_snp.R
#
# Writes process/AGE_SY/AGE_SY_4snpscan_no89.txt.gz
#   chr  pos  cM  sugar  sex  Wald_log10p  n_informative_founders
#
# The SNP scan is one row per SNP, so this is far larger than the haplotype
# equivalent -- the script prints the row count and file size, and names the scp
# command, so the size is known before anything is transferred.

suppressMessages(library(tidyverse))

DIR <- "process/AGE_SY"
OUT <- file.path(DIR, "AGE_SY_4snpscan_no89.txt.gz")

cells <- expand_grid(sugar = c("SY10", "SY20"), sex = c("F", "M")) %>%
  mutate(scan = paste0("AGE_", sugar, "_", sex, "_no89"),
         file = file.path(DIR, "Scans", scan, paste0(scan, ".snp_scan.txt")))

missing <- cells$file[!file.exists(cells$file)]
if (length(missing)) stop("missing:\n  ", paste(missing, collapse = "\n  "))

long <- cells %>% pmap_dfr(function(sugar, sex, scan, file) {
  d <- read.table(file, header = TRUE) %>% as_tibble()
  if (!"cM" %in% names(d)) d$cM <- NA_real_
  if (!"n_informative_founders" %in% names(d)) d$n_informative_founders <- NA_integer_
  d %>% transmute(chr, pos = as.integer(pos), cM, sugar, sex,
                  Wald_log10p, n_informative_founders)
})

write.table(long, gzfile(OUT), row.names = FALSE, quote = FALSE, sep = "\t")

cat("Wrote:", OUT, "\n")
cat(sprintf("  %s rows = %s SNPs x %d treatments\n",
            format(nrow(long), big.mark = ","),
            format(n_distinct(paste(long$chr, long$pos)), big.mark = ","),
            nrow(cells)))
cat(sprintf("  file size: %s\n",
            format(structure(file.size(OUT), class = "object_size"), units = "auto")))
long %>% group_by(sugar, sex) %>%
  summarise(snps = n(), max_wald = round(max(Wald_log10p, na.rm = TRUE), 1),
            median_informative = median(n_informative_founders, na.rm = TRUE),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nfetch with:\n")
cat("  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/",
    OUT, " process/AGE_SY/\n", sep = "")
