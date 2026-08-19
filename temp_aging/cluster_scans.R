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

suppressMessages({library(tidyverse)})

GATHER <- "process/AGE_2024/AGE_2024_vs_AGE_SY.txt.gz"
OUT    <- "temp_aging/FigureY_tree.png"

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

w <- g %>% inner_join(tiles, by = c("chr", "pos")) %>%
  mutate(id = paste0(diet, "_", reps)) %>%
  select(chr, pos, id, Wald_log10p) %>%
  pivot_wider(names_from = id, values_from = Wald_log10p) %>% drop_na()

m <- as.matrix(w[, names(NREP)])
cr <- cor(m, method = "spearman")
d  <- as.dist(1 - cr)
h  <- hclust(d, method = "average")

cat("Spearman correlation, ", nrow(m), " euchromatic tiles\n\n", sep = "")
print(round(cr, 2))

lab <- sprintf("%s  (n=%d)", h$labels, NREP[h$labels])

# replicate counts go on the tips: the axis is correlation, so a noisy n=2 scan
# sits far left for reasons that have nothing to do with food or cage.
png(OUT, width = 5.6, height = 3.6, units = "in", res = 300)
par(mar = c(4.2, 1, 1.5, 9), xpd = NA)
h2 <- h; h2$labels <- lab
plot(as.dendrogram(h2), horiz = TRUE, axes = FALSE,
     edgePar = list(lwd = 1.6, col = "grey25"))
at <- pretty(range(h$height), 5)
axis(1, at = at, labels = sprintf("%.2f", 1 - at), cex.axis = 0.8)
mtext("Spearman correlation", side = 1, line = 2.5, cex = 0.9)
dev.off()

cat("\nwrote", OUT, "\n\nmerge heights, as correlations:\n")
for (i in seq_len(nrow(h$merge))) {
  nm <- function(k) if (k < 0) h$labels[-k] else sprintf("[cluster %d]", k)
  cat(sprintf("  %-28s + %-28s  at r = %.2f\n",
              nm(h$merge[i,1]), nm(h$merge[i,2]), 1 - h$height[i]))
}
