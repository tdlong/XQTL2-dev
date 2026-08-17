# make_zoom_means.R — pull the founder frequencies around the seven zoom peaks.
#
# RUN THIS ON HPC3, from /dfs7/adl/tdlong/fly_pool/XQTL2-dev:
#   module load R/4.2.2
#   Rscript temp_aging/make_zoom_means.R
#
# Then fetch the single file it writes:
#   scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/process/AGE_SY/AGE_SY_zoom_means.txt.gz \
#       process/AGE_SY/
#
# The meansBySample.txt files are ~257 MB each and there are four of them. The
# zoom plots only need 1.2 Mb around each of seven peaks, which is well under 1%
# of that, so the subsetting happens on the cluster and one small file comes back.

suppressMessages(library(tidyverse))

SCANDIR <- "process/AGE_SY/Scans"
OUT     <- "process/AGE_SY/AGE_SY_zoom_means.txt.gz"
HALF    <- 0.6e6      # keep +/- this around each peak; plots use the middle 1 Mb

# Peak positions are the local Wald maxima nearest the positions Tony listed
# (X@10.1, 2L@10.6/12/15, 2R@14/16.5, 3L@9.2), taken from AGE_SY_4scan.txt.gz.
PEAKS <- tribble(
  ~locus,        ~chr,     ~pos,
  "chrX:10.09",  "chrX",   10090000,
  "chr2L:10.71", "chr2L",  10705000,
  "chr2L:11.95", "chr2L",  11945000,
  "chr2L:14.85", "chr2L",  14850000,
  "chr2R:13.88", "chr2R",  13880000,
  "chr2R:16.52", "chr2R",  16515000,
  "chr3L:9.31",  "chr3L",   9305000)

cells <- expand_grid(sugar = c("SY10", "SY20"), sex = c("F", "M")) %>%
  mutate(scan = paste0("AGE_", sugar, "_", sex),
         file = file.path(SCANDIR, scan, paste0(scan, ".meansBySample.txt")))

missing <- cells$file[!file.exists(cells$file)]
if (length(missing)) stop("missing:\n  ", paste(missing, collapse = "\n  "))

out <- pmap_dfr(cells, function(sugar, sex, scan, file) {
  d <- read.table(file, header = TRUE) %>% as_tibble()
  # one pass per peak, then bind: the file is long (chr, pos, TRT, REP, founder,
  # freq) so a single filter across all seven windows is simpler than a join
  pmap_dfr(PEAKS, function(locus, chr, pos) {
    d %>% filter(chr == !!chr, pos >= !!pos - HALF, pos <= !!pos + HALF) %>%
      mutate(locus = locus, peak_pos = !!pos)
  }) %>% mutate(sugar = sugar, sex = sex)
})

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.table(out, gzfile(OUT), row.names = FALSE, quote = FALSE, sep = "\t")

cat("Wrote:", OUT, "\n")
cat(sprintf("  %d rows; %d loci x %d treatments x %d founders x %d TRT x %d REP\n",
            nrow(out), n_distinct(out$locus), nrow(cells), n_distinct(out$founder),
            n_distinct(out$TRT), n_distinct(out$REP)))
out %>% count(locus, name = "rows") %>% as.data.frame() %>% print(row.names = FALSE)
