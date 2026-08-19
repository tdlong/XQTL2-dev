# cluster_scans.R — cluster the six food x cage Wald profiles.
#
# Six scans, females: SY10 / SY20 / lab food, each from the May 2023 and the
# November 2023 cage. The question is whether they group by food or by cage.
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript temp_aging/cluster_scans.R
#
# Distance is 1 - Spearman correlation over euchromatic non-overlapping windows,
# average linkage. Spearman rather than Pearson because the profiles have very
# different scales -- replicate counts run 2, 2, 4, 10, 10, 2, and the Wald
# statistic grows with them, so ranks are the comparable thing. Replicate counts
# are on the tip labels for the same reason: any tip with n=2 is noisy, which
# pushes it away from everything, cage or no cage.
#
# Dividing each scan by its replicate count, as FigureY does, changes NOTHING
# here: Spearman is rank-based and scaling by a positive constant leaves ranks
# alone. These are the raw scans' correlations either way.
#
# One tree per region. chr3R and the pericentromeric blocks are dominated by the
# heterochromatin rise, and chr3L by the major locus at ~9.3 Mb, so those are
# split out rather than left to drive everything.

suppressMessages({library(tidyverse)})

GATHER <- "process/AGE_2024/AGE_2024_vs_AGE_SY.txt.gz"
OUT    <- "temp_aging/FigureY_trees.png"

CHRS <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
HET <- tribble(~chr, ~eu_start, ~eu_end,
  "chrX", 2.5, 21.2, "chr2L", 0.5, 22.9, "chr2R", 1.3, 25.1,
  "chr3L", 0.7, 24.0, "chr3R", 4.5, 32.0)
NREP <- c(SY10_May = 2, SY20_May = 2, lab_May = 4,
          SY10_Nov = 10, SY20_Nov = 10, lab_Nov = 2)

g <- read.table(GATHER, header = TRUE, sep = "\t") %>% as_tibble() %>%
  filter(reps %in% c("May", "Nov"), chr %in% CHRS) %>%
  left_join(HET, by = "chr") %>%
  filter(pos/1e6 >= eu_start, pos/1e6 <= eu_end)

# every 15th window: 75 kb windows stepping 5 kb overlap fifteenfold
tiles <- g %>% distinct(chr, pos) %>% arrange(chr, pos) %>%
  group_by(chr) %>% slice(seq(1, n(), by = 15)) %>% ungroup()

base <- g %>% inner_join(tiles, by = c("chr", "pos")) %>%
  mutate(id = paste0(diet, "_", reps))

REGIONS <- list(
  "chr2L"        = function(d) d %>% filter(chr == "chr2L"),
  "chr3L"        = function(d) d %>% filter(chr == "chr3L"),
  "chr3L, major locus removed" =
    function(d) d %>% filter(chr == "chr3L", !(pos/1e6 > 8.5 & pos/1e6 < 10.5)),
  "whole euchromatic genome" = function(d) d)

fit <- function(sel) {
  w <- sel(base) %>% select(chr, pos, id, Wald_log10p) %>%
    pivot_wider(names_from = id, values_from = Wald_log10p) %>% drop_na()
  cr <- cor(as.matrix(w[, names(NREP)]), method = "spearman")
  list(cr = cr, n = nrow(w), h = hclust(as.dist(1 - cr), method = "average"))
}

res <- imap(REGIONS, function(f, nm) {
  r <- fit(f)
  cat("\n", nm, "  (", r$n, " tiles)\n\n", sep = "")
  print(round(r$cr, 2))
  r
})

png(OUT, width = 8.6, height = 3.4, units = "in", res = 300)
par(mfrow = c(1, 3), mar = c(4.2, 1, 2.2, 8.5), xpd = NA)
for (nm in c("chr2L", "chr3L", "chr3L, major locus removed")) {
  h <- res[[nm]]$h
  h2 <- h; h2$labels <- sprintf("%s  (n=%d)", h$labels, NREP[h$labels])
  plot(as.dendrogram(h2), horiz = TRUE, axes = FALSE,
       edgePar = list(lwd = 1.6, col = "grey25"))
  at <- pretty(range(h$height), 4)
  axis(1, at = at, labels = sprintf("%.2f", 1 - at), cex.axis = 0.9)
  mtext("Spearman correlation", side = 1, line = 2.6, cex = 0.7)
  mtext(nm, side = 3, line = 0.5, cex = 0.85, font = 2)
}
dev.off()
cat("\nwrote", OUT, "\n")

for (nm in names(res)) {
  h <- res[[nm]]$h
  cat("\n", nm, " -- merge heights as correlations:\n", sep = "")
  for (i in seq_len(nrow(h$merge))) {
    lbl <- function(k) if (k < 0) h$labels[-k] else sprintf("[%d]", k)
    cat(sprintf("  %-10s + %-10s  r = %.2f\n",
                lbl(h$merge[i,1]), lbl(h$merge[i,2]), 1 - h$height[i]))
  }
}
