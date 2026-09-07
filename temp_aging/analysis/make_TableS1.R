# make_TableS1.R — Supplementary Table 1: the well-separated peaks.
#
#   Rscript temp_aging/make_TableS1.R
#
# Writes temp_aging/TableS1.tsv, one row per peak, no commentary. peak_table.R
# does the peak finding and writes its working output to peak_table.txt; this
# turns that into the table the manuscript cites, adding the per-treatment
# heritability at each peak window from the scan.

suppressMessages(library(tidyverse))

PEAKS <- "temp_aging/peak_table.txt"   # paths are relative to the repo root
SCAN  <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
OUT   <- "temp_aging/TableS1.tsv"
for (f in c(PEAKS, SCAN)) if (!file.exists(f)) stop("missing ", f)

pk <- read.table(PEAKS, header = TRUE, sep = "\t", comment.char = "#") %>% as_tibble() %>%
  filter(grepl("^[AB]:", set))

sc <- read.table(SCAN, header = TRUE, sep = "\t") %>% as_tibble() %>%
  mutate(sex = ifelse(sex %in% c("FALSE", "F"), "F", "M"),
         trt = paste0(sugar, "_", sex))

# h2 at each peak window, per treatment
h2 <- pk %>% transmute(chr, pos = round(Mb * 1e6)) %>%
  left_join(sc %>% select(chr, pos, trt, H2), by = c("chr", "pos")) %>%
  pivot_wider(names_from = trt, values_from = H2, names_prefix = "h2_")

ARMS <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

body <- pk %>%
  mutate(pos = round(Mb * 1e6)) %>%
  left_join(h2, by = c("chr", "pos")) %>%
  mutate(.A = grepl("^A", set), .arm = factor(chr, levels = ARMS)) %>%
  arrange(desc(.A), .arm, Mb) %>%
  transmute(
    .A,
    Chromosome        = chr,
    `Position (Mb)`   = sprintf("%.2f", Mb),
    `Position (cM)`   = sprintf("%.1f", cM),
    `Interval (Mb)`   = sprintf("%.2f-%.2f", int_from_Mb, int_to_Mb),
    `Width (kb)`      = as.character(int_kb),
    `Width (cM)`      = sprintf("%.2f", int_cM),
    `Peak treatment`  = peak_trait,
    `-log10P SY10 F`  = sprintf("%.1f", SY10_F),
    `-log10P SY20 F`  = sprintf("%.1f", SY20_F),
    `-log10P SY10 M`  = sprintf("%.1f", SY10_M),
    `-log10P SY20 M`  = sprintf("%.1f", SY20_M),
    `h2 SY10 F`       = sprintf("%.2f", h2_SY10_F),
    `h2 SY20 F`       = sprintf("%.2f", h2_SY20_F),
    `h2 SY10 M`       = sprintf("%.2f", h2_SY10_M),
    `h2 SY20 M`       = sprintf("%.2f", h2_SY20_M))

# a subtitle row: the label in the first column, the rest blank
subtitle <- function(txt, template) {
  r <- template[1, , drop = FALSE]
  r[1, ] <- ""
  r[1, 1] <- txt
  r
}
cols <- setdiff(names(body), ".A")
A <- body %>% filter(.A)  %>% select(all_of(cols))
B <- body %>% filter(!.A) %>% select(all_of(cols))

tab <- bind_rows(
  subtitle(sprintf("Peaks exceeding -log10 P of 15 (n = %d)", nrow(A)), A), A,
  subtitle(sprintf("Peaks between -log10 P of 7.5 and 15 (n = %d)", nrow(B)), B), B)

write.table(tab, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", OUT, "--", nrow(A), "peaks above 15,", nrow(B), "between 7.5 and 15\n\n")
print(as.data.frame(tab), row.names = FALSE)
