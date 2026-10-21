# drift_tiles.R — per-replicate founder frequencies, control and selected pools, on a
# 0.25 cM grid genome-wide. Input to a LOCAL LD score per tile (ld_drift.R
# measures the same thing in seven 1.2 Mb regions and finds it varies by region).
#
#   sbatch temp_aging/analysis/drift_tiles.sh        # HPC3, reads meansBySample
#
# Writes temp_aging/numbers/drift_tiles.txt.gz: chr, cM, pos, sex, sugar, pool
# (C or Z), REP, founder, freq. Controls are one pool per replicate per sex; they are
# written once, under sugar = "both". About 700k rows.

suppressMessages(library(tidyverse))

SCANDIR  <- "process/AGE_SY/Scans"
DESDIR   <- "helpfiles/AGE_SY"
CMSRC    <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
OUT      <- "temp_aging/numbers/drift_tiles.txt.gz"
DROP_REP <- c(8, 9)
STEP     <- 0.25

cm <- read.table(CMSRC, header = TRUE, sep = "\t") %>% as_tibble() %>% distinct(chr, pos, cM)
# one window per 0.25 cM: the window nearest each grid point
grid <- cm %>% mutate(g = round(cM / STEP) * STEP) %>% group_by(chr, g) %>%
  slice_min(abs(cM - g), n = 1, with_ties = FALSE) %>% ungroup() %>% select(chr, pos, cM = g)
cat("grid points:", nrow(grid), "\n")

cells <- expand_grid(sugar = c("SY10", "SY20"), sex = c("F", "M")) %>%
  mutate(scan = paste0("AGE_", sugar, "_", sex),
         file = file.path(SCANDIR, scan, paste0(scan, ".meansBySample.txt")))
miss <- cells$file[!file.exists(cells$file)]
if (length(miss)) stop("missing:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)

out <- list()
for (r in seq_len(nrow(cells))) {
  cl <- cells[r, ]; cat("==", cl$scan, "\n")
  d <- read.table(cl$file, header = TRUE) %>% as_tibble() %>%
    filter(!REP %in% DROP_REP, TRT %in% c("C", "Z")) %>%
    inner_join(grid, by = c("chr", "pos")) %>%
    transmute(chr, cM, pos, sex = cl$sex, sugar = if_else(TRT == "C", "both", cl$sugar),
              pool = TRT, REP, founder, freq)
  if (cl$sugar == "SY20") d <- d %>% filter(pool == "Z")   # controls already written under SY10
  out[[r]] <- d
}
out <- bind_rows(out)
cat("rows:", nrow(out), "\n")
con <- gzfile(OUT, "w"); write.table(out, con, row.names = FALSE, quote = FALSE, sep = "\t"); close(con)
cat("Wrote:", OUT, "\n")
