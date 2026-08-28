# scan_resolution.R — the genome in 1 cM tiles: how much responds, how wide peaks are.
#
#   Rscript temp_aging/scan_resolution.R
#
# DATA: the 10-replicate AGE_SY dataset -- replicates 8 and 9 dropped, those
# being the May 2023 cage (helpfiles/AGE_2024/population_assignment.txt).
#
# TILES ARE 1 cM. Mapping resolution is genetic, so a genetic tile is the unit
# in which the genome divides into roughly independent pieces. A physical tile
# is not: 1 Mb is a fraction of a cM near a centromere and several cM in the
# middle of an arm, so a physical tile counts low-recombination sequence many
# times over and high-recombination sequence too few.
#
# The euchromatic map is ~267 cM, so this is ~267 tiles.

suppressMessages(library(tidyverse))

TILE <- 1                       # cM
d <- read.table("process/AGE_SY/AGE_SY_4scan_no89.txt.gz", header = TRUE,
                sep = "\t") %>% as_tibble()
HET <- read.table("pipeline/helpfiles/het_bounds.txt", header = TRUE,
                  comment.char = "#") %>% as_tibble()

eu <- d %>% left_join(HET, by = "chr") %>%
  filter(pos/1e6 >= eu_start, pos/1e6 <= eu_end)

# one row per 1 cM tile per arm, carrying the strongest window in it
tl <- eu %>%
  group_by(chr, pos, cM) %>% summarise(wald = max(Wald_log10p), .groups = "drop") %>%
  mutate(tile = floor(cM / TILE)) %>%
  group_by(chr, tile) %>%
  summarise(cM_lo = min(cM), cM_hi = max(cM), Mb = (max(pos) - min(pos))/1e6,
            wald = max(wald), n_win = n(), .groups = "drop")

cat(sprintf("%d tiles of %g cM over %d cM of euchromatin (%.0f Mb)\n\n",
            nrow(tl), TILE, round(sum(tl$cM_hi - tl$cM_lo)), sum(tl$Mb)))

cat("fraction of tiles whose strongest window exceeds a threshold\n")
cat("(max over the four treatments)\n\n")
map_dfr(c(5, 10, 15, 20), function(th)
  tibble(threshold = th, n = sum(tl$wald > th),
         `% of tiles` = round(100 * mean(tl$wald > th), 1))) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\ntile size in Mb -- the reason a physical tile is the wrong unit\n\n")
tl %>% summarise(`median Mb` = round(median(Mb), 2),
                    `min Mb` = round(min(Mb), 2),
                    `max Mb` = round(max(Mb), 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nby arm:\n\n")
tl %>% group_by(chr) %>%
  summarise(n = n(), `% > 5` = round(100*mean(wald > 5), 0),
            `% > 15` = round(100*mean(wald > 15), 0),
            `median Mb/tile` = round(median(Mb), 2), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\ntiles above 15, in order:\n\n")
tl %>% filter(wald > 15) %>% arrange(chr, tile) %>%
  transmute(chr, cM = paste0(round(cM_lo,1), "-", round(cM_hi,1)),
            Mb = round(Mb, 2), wald = round(wald, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)
