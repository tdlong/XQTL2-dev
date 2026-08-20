# gather.R — collapse the twelve no-8-9 scans into the two files the figures read.
#
# RUN ON HPC3 from the repo root, once run_scans.sh has finished:
#   module load R/4.2.2
#   Rscript scripts_oneoffs/AGE_SY/nov_only/gather.R
#
# Writes, in the same format as their 12-replicate counterparts so the figure
# scripts read them unchanged:
#
#   process/AGE_SY/AGE_SY_4scan_no89.txt.gz
#       chr pos cM sugar sex Wald_log10p Cutl_H2 Falc_H2 <bias cols>
#       -- panels A and B of Figure 1, and the Wald panels of Figure 2
#
#   process/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz
#       chr pos sex sugar half Cutl_H2 Falc_H2 Wald_log10p <bias cols>
#       -- the split-half error term, hence panel C of Figure 1
#
# Fetch both:
#   scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/process/AGE_SY/AGE_SY_4scan_no89.txt.gz process/AGE_SY/
#   scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/process/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz process/AGE_SY_splithalf/

suppressMessages(library(tidyverse))

FULL_DIR <- "process/AGE_SY"
HALF_DIR <- "process/AGE_SY_splithalf"
OUT_FULL <- file.path(FULL_DIR, "AGE_SY_4scan_no89.txt.gz")
OUT_HALF <- file.path(HALF_DIR, "AGE_SY_splithalf_H2_no89.txt.gz")

# bias columns exist only in scans made after XQTL2 #34; cM comes from the
# pipeline's add_genetic(). Carry whichever are present, NA otherwise.
OPTIONAL <- c("Cutl_H2_bias", "Falc_H2_bias", "Cutl_clamp_frac", "cM")
fill <- function(d) { for (v in OPTIONAL) if (!v %in% names(d)) d[[v]] <- NA_real_; d }

read_scan <- function(dir, scan) {
  f <- file.path(dir, "Scans", scan, paste0(scan, ".scan.txt"))
  if (!file.exists(f)) { cat("  MISSING:", f, "\n"); return(NULL) }
  read.table(f, header = TRUE) %>% as_tibble() %>% fill()
}

cells <- expand_grid(sugar = c("SY10", "SY20"), sex = c("F", "M")) %>%
  mutate(scan = paste0("AGE_", sugar, "_", sex))

# ── the four 10-replicate scans ──────────────────────────────────────────────
cat("10-replicate scans:\n")
full <- cells %>% pmap_dfr(function(sugar, sex, scan) {
  d <- read_scan(FULL_DIR, paste0(scan, "_no89")); if (is.null(d)) return(NULL)
  d %>% transmute(chr, pos = as.integer(pos), cM, sugar, sex,
                  Wald_log10p, Cutl_H2, Falc_H2,
                  Cutl_H2_bias, Falc_H2_bias, Cutl_clamp_frac)
})
if (!nrow(full)) stop("none of the four 10-replicate scans were found")

# ── the eight halves ─────────────────────────────────────────────────────────
cat("5-replicate halves:\n")
half <- expand_grid(sugar = c("SY10", "SY20"), sex = c("F", "M"),
                    half = c("odd", "even")) %>%
  pmap_dfr(function(sugar, sex, half) {
    scan <- paste0("AGE_", sugar, "_", sex, "_no89_", half)
    d <- read_scan(HALF_DIR, scan); if (is.null(d)) return(NULL)
    d %>% transmute(chr, pos = as.integer(pos), sex, sugar, half,
                    Cutl_H2, Falc_H2, Wald_log10p,
                    Cutl_H2_bias, Falc_H2_bias, Cutl_clamp_frac)
  })
if (!nrow(half)) stop("none of the eight halves were found")

dir.create(FULL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(HALF_DIR, showWarnings = FALSE, recursive = TRUE)
write.table(full, gzfile(OUT_FULL), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(half, gzfile(OUT_HALF), row.names = FALSE, quote = FALSE, sep = "\t")

cat("\nWrote:\n  ", OUT_FULL, "\n  ", OUT_HALF, "\n\n", sep = "")
full %>% group_by(sugar, sex) %>%
  summarise(windows = n(), max_wald = round(max(Wald_log10p, na.rm = TRUE), 1),
            median_h2 = round(median(Cutl_H2, na.rm = TRUE), 3), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nhalves:\n")
half %>% count(sugar, sex, half, name = "windows") %>%
  as.data.frame() %>% print(row.names = FALSE)
