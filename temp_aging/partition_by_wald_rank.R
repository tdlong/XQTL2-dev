# partition_by_wald_rank.R — how much of the heritability sits in the strongest
# part of the genome, and how much of that is sex + diet.
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript temp_aging/partition_by_wald_rank.R
#
# Blocks are 2 cM, NON-OVERLAPPING, euchromatin only. Each of those matters:
#
#   cM not Mb   h2 is smeared over ~9 Mb in the low-recombination regions
#               flanking the centromeres and ~1 Mb in euchromatin, so a fixed
#               physical block counts one region as ten. A 2 cM block here runs
#               0.38 to 6.38 Mb depending on where it sits.
#   euchromatin heterochromatin is 7% of the sequence but carries 20% of the h2
#               -- a threefold enrichment -- and reads ~0% sex+diet. Dropping it
#               raises the estimate.
#   no overlap  the scan uses 75 kb windows stepping 5 kb, so windows overlap
#               15-fold. Every 15th window tiles the genome exactly once;
#               summing raw windows counts each piece of genome 15 times.
#
# The h2 column is scaled so the whole genome sums to 50%, a nominal heritability
# for longevity. That scaling is a factor of ~6 and assumes the overcounting is
# uniform, so the ABSOLUTE columns are indicative only. The RATIO does not depend
# on the scaling at all.

suppressMessages({library(tidyverse); library(splines)})

IN     <- "process/AGE_SY_splithalf/H2_varcomp_by_window.txt.gz"
FLYMAP <- "pipeline/helpfiles/flymap.r6.txt"   # symlink to the XQTL2 install
CM     <- 2      # block width, cM
H2_TOT <- 50     # assumed trait heritability, for the absolute columns

# dm6 euchromatin boundaries (Huynh et al. 2023 PLoS Genet 19:e1010439, Table S2)
HET <- tribble(~chr, ~eu_start, ~eu_end,
  "chr2L", 0.5, 22.9, "chr2R", 1.3, 25.1, "chr3L", 0.7, 24.0, "chr3R", 4.5, 32.0)

# same conversion as add_genetic() in the pipeline's scan_functions.R
add_genetic <- function(df) {
  fm <- read.table(FLYMAP, header = FALSE); colnames(fm) <- c("chr","pos","cM")
  df$cM <- NA_real_
  for (c_ in unique(df$chr)) {
    x <- fm %>% filter(chr == c_)
    o <- ksmooth(x$pos, x$cM, kernel = "normal", bandwidth = 3e6)
    df$cM[df$chr == c_] <- splinefun(o$x, o$y)(df$pos[df$chr == c_])
  }
  df
}

v <- read.table(IN, header = TRUE, sep = "\t") %>% as_tibble() %>%
  filter(chr != "chrX") %>%                       # X excluded: see SUMMARY.md
  mutate(other = sex + diet + sex.diet) %>%
  add_genetic() %>%
  left_join(HET, by = "chr") %>%
  filter(pos/1e6 >= eu_start, pos/1e6 <= eu_end)

# every 15th window tiles the genome once (75 kb windows, 5 kb steps)
til <- v %>% arrange(chr, pos) %>% group_by(chr) %>% slice(seq(1, n(), by = 15)) %>% ungroup()

blk <- til %>%
  mutate(b = paste0(chr, "_", floor(cM / CM))) %>%
  group_by(chr, b) %>%
  summarise(Mb = n() * 0.075, cM_lo = min(cM), pos = min(pos),
            H2 = sum(H2), main = sum(main), other = sum(other),
            wald = max(wald_max), tiles = n(), .groups = "drop") %>%
  filter(tiles >= 4) %>%
  mutate(tot = main + other, fo = ifelse(tot > 0, other/tot, NA),
         h2_other = H2 * fo, h2_main = H2 * (1 - fo),
         grp = paste0(chr, "_", floor(cM_lo / 8)))   # 8 cM groups, bootstrap unit

K <- H2_TOT / sum(blk$H2)
cat(sprintf("euchromatin: %d blocks of %g cM covering %.0f Mb\n", nrow(blk), CM, sum(blk$Mb)))
cat(sprintf("block size in Mb: median %.2f, range %.2f-%.2f\n",
            median(blk$Mb), min(blk$Mb), max(blk$Mb)))
cat(sprintf("summed h2 %.0f points, scaled x%.3f so the genome totals %g%%\n\n",
            sum(blk$H2), K, H2_TOT))

set.seed(1)
res <- map_dfr(c(5, 10, 20, 30, 50, 100), function(pc) {
  s  <- blk %>% filter(wald >= quantile(blk$wald, max(0, 1 - pc/100)))
  g  <- unique(s$grp); ix <- split(seq_len(nrow(s)), s$grp)
  bs <- replicate(4000, { x <- s[unlist(ix[sample(g, length(g), TRUE)]), ]
                          100 * sum(x$h2_other, na.rm = TRUE) / sum(x$H2) })
  tibble(`top % genome` = pc, blocks = nrow(s), Mb = round(sum(s$Mb), 1),
         `h2 (of 50)` = round(K * sum(s$H2), 1),
         main = round(K * sum(s$h2_main, na.rm = TRUE), 1),
         `sex+diet` = round(K * sum(s$h2_other, na.rm = TRUE), 2),
         `% sex+diet` = round(100 * sum(s$h2_other, na.rm = TRUE) / sum(s$H2), 1),
         CI = sprintf("[%.1f, %.1f]", quantile(bs, .025, na.rm = TRUE),
                                      quantile(bs, .975, na.rm = TRUE)))
})
print(as.data.frame(res), row.names = FALSE)

cat("\nThe 100% row is not a conservative choice -- it is broken. In the bottom\n")
cat("half of the genome the non-main terms sit BELOW the replicate error, so they\n")
cat("sum negative and drag the numerator through zero.\n")
