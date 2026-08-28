# peak_table.R — the well-separated peaks, with the trait each belongs to.
#
#   Rscript temp_aging/peak_table.R
#
# DATA: the 10-replicate AGE_SY dataset -- replicates 8 and 9 dropped, those
# being the May 2023 cage (helpfiles/AGE_2024/population_assignment.txt).
#
# DEFINITION, since the text has to state it:
#   - work with the maximum -log10 P over the four treatments at each window
#   - take the tallest window in the euchromatin; that is a peak
#   - exclude EXCL cM either side of it on that arm; repeat while the running
#     maximum still exceeds THRESH
#   - the 2-LOD interval is read off the scan itself: from the peak window, walk
#     left and right while the value stays within 2 of the peak. Those are
#     PHYSICAL boundaries; the cM width is those same boundaries converted.
#
# EXCL is a choice, not a property of the data, and it sets how many peaks there
# are. It is written into the output so the number can never be quoted without it.
#
# Writes temp_aging/peak_table.txt (and prints it).

suppressMessages(library(tidyverse))

SCAN   <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"
OUT    <- "temp_aging/peak_table.txt"
THRESH <- 15      # -log10 P a peak must exceed
EXCL   <- 5       # cM excluded either side of a peak before looking for the next
DROP   <- 2       # -log10 P drop defining the support interval

d <- read.table(SCAN, header = TRUE, sep = "\t") %>% as_tibble()
HET <- read.table("pipeline/helpfiles/het_bounds.txt", header = TRUE,
                  comment.char = "#") %>% as_tibble()
eu <- d %>% left_join(HET, by = "chr") %>%
  filter(pos/1e6 >= eu_start, pos/1e6 <= eu_end)

# per window: the strongest treatment, and which one it is
w <- eu %>% group_by(chr, pos, cM) %>%
  slice_max(Wald_log10p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(chr, pos, cM, wald = Wald_log10p,
            trait = paste0(sugar, " ", ifelse(sex == "F", "female", "male")))

# top-down with an exclusion radius
pk <- list(); x <- w
while (nrow(x) && max(x$wald) > THRESH) {
  p <- x[which.max(x$wald), ]; pk[[length(pk) + 1]] <- p
  x <- x %>% filter(!(chr == p$chr & abs(cM - p$cM) <= EXCL))
}
pk <- bind_rows(pk)

# support interval, walked in PHYSICAL order, then converted to cM
tab <- pk %>% pmap_dfr(function(chr, pos, cM, wald, trait) {
  y <- w %>% filter(chr == !!chr) %>% arrange(pos)
  k <- which(y$pos == pos); lim <- wald - DROP
  lo <- k; while (lo > 1        && y$wald[lo - 1] > lim) lo <- lo - 1
  up <- k; while (up < nrow(y)  && y$wald[up + 1] > lim) up <- up + 1
  # every treatment's value at the peak window, so the table shows whether the
  # peak is specific to one trait or shared
  at <- eu %>% filter(chr == !!chr, pos == !!pos) %>%
    transmute(t = paste0(sugar, "_", ifelse(sex == "F", "F", "M")),
              v = round(Wald_log10p, 1)) %>% deframe()
  tibble(chr, Mb = round(pos/1e6, 2), cM = round(cM, 1),
         peak_trait = trait, peak_wald = round(wald, 1),
         SY10_F = at[["SY10_F"]], SY20_F = at[["SY20_F"]],
         SY10_M = at[["SY10_M"]], SY20_M = at[["SY20_M"]],
         int_cM = round(y$cM[up] - y$cM[lo], 2),
         int_kb = round((y$pos[up] - y$pos[lo]) / 1e3),
         int_from_Mb = round(y$pos[lo]/1e6, 2), int_to_Mb = round(y$pos[up]/1e6, 2))
}) %>% arrange(desc(peak_wald))

hdr <- sprintf(paste0(
  "# Well-separated peaks, AGE_SY 10-replicate dataset (replicates 8, 9 dropped)\n",
  "# peak       = tallest window, taken top-down over the max across the four traits\n",
  "# separated  = at least %g cM from any stronger peak on the same arm\n",
  "# threshold  = -log10 P > %g\n",
  "# interval   = walk left/right from the peak while within %g of it; physical,\n",
  "#              then converted to cM\n",
  "# peak_trait = which treatment is highest AT the peak window\n",
  "# SY10_F ... = every treatment's -log10 P at that same window\n",
  "# source     = %s\n#\n"), EXCL, THRESH, DROP, SCAN)
con <- file(OUT, "w"); writeLines(hdr, con)
suppressWarnings(write.table(tab, con, sep = "\t", quote = FALSE, row.names = FALSE))
close(con)

cat(hdr)
as.data.frame(tab) %>% print(row.names = FALSE)
cat(sprintf("\n%d peaks. interval: median %.2f cM, range %.2f-%.2f; %.0f kb to %.2f Mb\n",
            nrow(tab), median(tab$int_cM), min(tab$int_cM), max(tab$int_cM),
            min(tab$int_kb), max(tab$int_kb)/1e3))
cat("wrote", OUT, "\n")
