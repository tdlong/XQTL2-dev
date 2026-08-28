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

cat("h2 at a given Wald, per trait (spline of h2 on log Wald, over windows)\n")
cat("and the fraction of 1 cM TILES whose strongest window exceeds it\n\n")
res <- map_dfr(sort(unique(w$trt)), function(t) {
  z <- w %>% filter(trt == t) %>% mutate(lw = log(Wald_log10p))
  m <- gam(Cutl_H2 ~ s(lw), data = z)
  p <- as.numeric(predict(m, newdata = tibble(lw = log(c(2, 5, 10)))))
  h5 <- p[2]
  tl <- z %>% group_by(chr, tile) %>%
    summarise(h2 = max(Cutl_H2), wald = max(Wald_log10p), .groups = "drop")
  tibble(trait = t,
         `h2 @W2` = round(p[1], 2), `h2 @W5` = round(h5, 2), `h2 @W10` = round(p[3], 2),
         `% tiles h2 > that` = round(100 * mean(tl$h2 > h5), 1),
         `% tiles Wald > 5`  = round(100 * mean(tl$wald > 5), 1))
})
as.data.frame(res) %>% print(row.names = FALSE)
cat(sprintf("\n  h2 at Wald 5 spans %.2f to %.2f across the four traits\n",
            min(res$`h2 @W5`), max(res$`h2 @W5`)))
cat(sprintf("  tiles exceeding it: %.0f%% to %.0f%%\n",
            min(res$`% tiles h2 > that`), max(res$`% tiles h2 > that`)))

cat("\nfloor: h2 at tiles whose own Wald stays below 2\n\n")
w %>% group_by(chr, tile, trt) %>%
  summarise(h2 = max(Cutl_H2), wald = max(Wald_log10p), .groups = "drop") %>%
  filter(wald < 2) %>% group_by(trt) %>%
  summarise(tiles = n(), `median h2` = round(median(h2), 2), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nh2 at the peaks in Table S1, per trait\n\n")
PK <- tribble(~chr, ~Mb,
  "chr3L", 9.30, "chr3R", 8.63, "chr3L", 21.62, "chr2L", 11.95,
  "chr2L", 14.85, "chr3L", 7.43, "chr2R", 13.88, "chr2R", 24.08)
PK %>% pmap_dfr(function(chr, Mb) {
  w %>% filter(chr == !!chr, abs(pos/1e6 - Mb) < 0.03) %>%
    group_by(trt) %>% summarise(h2 = max(Cutl_H2), .groups = "drop") %>%
    pivot_wider(names_from = trt, values_from = h2) %>%
    mutate(peak = paste0(chr, " ", Mb), .before = 1)
}) %>% mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\ndistribution of tile h2 (max over the four traits)\n\n")
tl <- w %>% group_by(chr, tile) %>% summarise(h2 = max(Cutl_H2), .groups = "drop")
print(round(quantile(tl$h2, c(.1,.25,.5,.75,.9,1)), 2))
