# h2_threshold.R — what h2 corresponds to a significant Wald, and how much of
# the genome exceeds it.
#
#   Rscript temp_aging/h2_threshold.R
#
# DATA: the 10-replicate AGE_SY dataset -- replicates 8 and 9 dropped, those
# being the May 2023 cage (helpfiles/AGE_2024/population_assignment.txt).
#
# Within a trait, Wald and h2 are driven by the same founder frequency shifts and
# trace one curve, so a threshold on the Wald scale transfers onto the h2 scale.
# The curve is fitted per trait over WINDOWS -- it is a relationship between two
# quantities, not a statement about how much genome does anything.
#
# The genome fractions ARE tiled: 1 cM tiles cut on the map, keeping any tile not
# wholly heterochromatic, matching scan_resolution.R. A tile's value is the
# largest of its windows, as there.
#
# h2 is raw and uncorrected, as plotted in Fig 1b.

suppressMessages({library(tidyverse); library(mgcv)})

d <- read.table("process/AGE_SY/AGE_SY_4scan_no89.txt.gz", header = TRUE,
                sep = "\t") %>% as_tibble()
HET <- read.table("pipeline/helpfiles/het_bounds.txt", header = TRUE,
                  comment.char = "#") %>% as_tibble()

w <- d %>% left_join(HET, by = "chr") %>% group_by(chr) %>%
  mutate(tile = floor((cM - min(cM)) / 1),
         is_eu = pos/1e6 >= eu_start & pos/1e6 <= eu_end) %>% ungroup() %>%
  mutate(trt = paste0(sugar, "_", sex))
keep <- w %>% group_by(chr, tile) %>% summarise(f = mean(is_eu), .groups = "drop") %>%
  filter(f > 0) %>% select(chr, tile)
w <- w %>% semi_join(keep, by = c("chr", "tile")) %>% filter(is_eu)

# The floor and the cutoffs, both on 1 cM tiles. The floor is the flat average
# over tiles whose Wald never reaches 2; the cutoffs are a smooth fit of tile h2
# on log(tile Wald) read off at 7.5 and 15. Same object, estimated the same way,
# so the three numbers are comparable.
#
# This is NOT a conversion. h2 and Wald are correlated but not monotone -- within
# a trait Spearman is 0.69-0.79, and among windows at Wald 7-8 the h2 runs 0.78
# to 1.35 between the 10th and 90th percentiles. These are averages through wide
# scatter: the h2 a tile at that significance carries on average, not the h2 that
# significance implies.
# A tile is read at ONE step -- the one with the largest Wald -- as in
# significant_regions.R and scan_resolution.R. Taking max(H2) and max(Wald)
# separately would let them come from different steps. It moves the numbers by
# ~0.003 percentage points, but the tile rule has to mean one thing.
tl <- w %>% mutate(trt = paste0(sugar, "_", sex)) %>%
  group_by(trt, chr, tile) %>%
  slice_max(Wald_log10p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(trt, chr, tile, h2 = H2, wald = Wald_log10p)

cat("h2 on 1 cM tiles: the floor, and the value at each cutoff\n\n")
res <- map_dfr(sort(unique(tl$trt)), function(t) {
  z <- tl %>% filter(trt == t) %>% mutate(lw = log(pmax(wald, 0.01)))
  m <- gam(h2 ~ s(lw), data = z)
  p <- as.numeric(predict(m, newdata = tibble(lw = log(c(7.5, 15)))))
  tibble(trait = t,
         `floor (Wald<2)` = round(mean(z$h2[z$wald < 2]), 2),
         `n tiles` = sum(z$wald < 2),
         `at 7.5` = round(p[1], 2), `at 15` = round(p[2], 2),
         `% tiles > 7.5` = round(100 * mean(z$wald > 7.5), 1),
         `% tiles > 15`  = round(100 * mean(z$wald > 15), 1))
})
as.data.frame(res) %>% print(row.names = FALSE)
cat(sprintf("\n  averaged over the four traits:  floor %.2f   at 7.5 %.2f   at 15 %.2f\n",
            mean(res$`floor (Wald<2)`), mean(res$`at 7.5`), mean(res$`at 15`)))
cat(sprintf("  ranges: %.2f-%.2f, %.2f-%.2f, %.2f-%.2f\n",
            min(res$`floor (Wald<2)`), max(res$`floor (Wald<2)`),
            min(res$`at 7.5`), max(res$`at 7.5`),
            min(res$`at 15`), max(res$`at 15`)))

cat("\nh2 at the peaks in Table S1, per trait\n\n")
PK <- tribble(~chr, ~Mb,
  "chr3L", 9.30, "chr3R", 8.63, "chr3L", 21.62, "chr2L", 11.95,
  "chr2L", 14.85, "chr3L", 7.43, "chr2R", 13.88, "chr2R", 24.08)
PK %>% pmap_dfr(function(chr, Mb) {
  w %>% filter(chr == !!chr, abs(pos/1e6 - Mb) < 0.03) %>%
    group_by(trt) %>% summarise(h2 = max(H2), .groups = "drop") %>%
    pivot_wider(names_from = trt, values_from = h2) %>%
    mutate(peak = paste0(chr, " ", Mb), .before = 1)
}) %>% mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\ndistribution of tile h2 (max over the four traits)\n\n")
tl <- w %>% group_by(chr, tile) %>% summarise(h2 = max(H2), .groups = "drop")
print(round(quantile(tl$h2, c(.1,.25,.5,.75,.9,1)), 2))
