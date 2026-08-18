# compare_depth.R — read depth per sample, pilot against AGE_SY.
#
# RUN ON HPC3 from the repo root:
#   module load R/4.2.2
#   Rscript scripts_oneoffs/AGE_2024/compare_depth.R
#
# The question is whether the Aug 2024 pilot is noisier than AGE_SY because it
# was sequenced shallower. Its h2 sits ~0.25 higher genome-wide including where
# nothing is happening, and it calls twice as many windows significant, both of
# which are what thinner coverage would do.
#
# Depth is measured where it matters: at the catalog SNPs actually used, from
# RefAlt.<chr>.txt, as REF+ALT summed per sample. chr2L only -- one chromosome is
# plenty for a per-sample mean and the files are large.

suppressMessages(library(tidyverse))

CHR <- "chr2L"
SETS <- tribble(~label,     ~file,
  "AGE_2024 (pilot)",  file.path("process/AGE_2024", "Calls", paste0("RefAlt.", CHR, ".txt")),
  "AGE_SY",            file.path("process/AGE_SY",   "Calls", paste0("RefAlt.", CHR, ".txt")))

depth_of <- function(file) {
  if (!file.exists(file)) { warning("missing ", file); return(NULL) }
  d <- read.table(file, header = TRUE, nrows = 2e5, check.names = FALSE)
  ref <- grep("^REF_", names(d), value = TRUE)
  tibble(sample = sub("^REF_", "", ref)) %>%
    mutate(mean_depth = map_dbl(sample, function(s)
      mean(d[[paste0("REF_", s)]] + d[[paste0("ALT_", s)]], na.rm = TRUE)))
}

res <- SETS %>% pmap_dfr(function(label, file) {
  x <- depth_of(file); if (is.null(x)) return(NULL); x %>% mutate(set = label)
})
if (!nrow(res)) stop("no RefAlt tables found")

# founders are in both and are not the comparison
res <- res %>% filter(!sample %in% c("B1","B2","B3","B4","B5","B6","B7","AB8"))

cat("mean read depth at", CHR, "catalog SNPs (first 200k rows)\n\n")
res %>% group_by(set) %>%
  summarise(samples = n(),
            median = round(median(mean_depth), 1),
            min = round(min(mean_depth), 1),
            max = round(max(mean_depth), 1), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nper sample, pilot:\n")
res %>% filter(set == "AGE_2024 (pilot)") %>% arrange(sample) %>%
  transmute(sample, depth = round(mean_depth, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nAGE_SY, the 6 female cages the pilot is compared against:\n")
res %>% filter(set == "AGE_SY", grepl("_R[1-6]_F$", sample)) %>% arrange(sample) %>%
  transmute(sample, depth = round(mean_depth, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)
