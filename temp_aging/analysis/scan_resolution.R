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
# TILING IS ON THE MAP, NOT ON THE DATA. Each arm is cut at 0-1, 1-2, 2-3 cM and
# so on, from the arm's own start. A tile is kept if ANY part of it lies inside
# the euchromatin boundaries (pipeline/helpfiles/het_bounds.txt); a tile lying
# wholly in heterochromatin is dropped. Tiles that straddle a boundary count --
# they are partly euchromatic, and the cut is on the map, not on where the
# heterochromatin happens to start.

suppressMessages(library(tidyverse))

TILE <- 1                       # cM
d <- read.table("process/AGE_SY/AGE_SY_4scan_no89.txt.gz", header = TRUE,
                sep = "\t") %>% as_tibble()
HET <- read.table("pipeline/helpfiles/het_bounds.txt", header = TRUE,
                  comment.char = "#") %>% as_tibble()

# cut every arm into 1 cM tiles from its own start, then mark each window with
# which tile it falls in and whether it is euchromatic
w <- d %>% left_join(HET, by = "chr") %>%
  group_by(chr) %>%
  mutate(tile = floor((cM - min(cM)) / TILE),
         is_eu = pos/1e6 >= eu_start & pos/1e6 <= eu_end) %>%
  ungroup()

tl <- w %>% group_by(chr, tile, pos, cM, is_eu) %>%
  summarise(wald = max(Wald_log10p), .groups = "drop") %>%
  group_by(chr, tile) %>%
  summarise(cM_lo = min(cM), cM_hi = max(cM), Mb = (max(pos) - min(pos))/1e6,
            wald = max(wald), n_win = n(),
            f_eu = mean(is_eu), .groups = "drop")

n_all <- nrow(tl)
n_het <- sum(tl$f_eu == 0)
n_str <- sum(tl$f_eu > 0 & tl$f_eu < 1)
tl <- tl %>% filter(f_eu > 0)          # keep anything touching euchromatin
cat(sprintf("%d tiles of %g cM across the whole map; %d lie wholly in\n",
            n_all, TILE, n_het))
cat(sprintf("heterochromatin and are dropped, %d straddle a boundary and are kept.\n\n",
            n_str))

cat(sprintf("%d tiles kept, spanning %d cM and %.0f Mb\n\n",
            nrow(tl), round(sum(tl$cM_hi - tl$cM_lo)), sum(tl$Mb)))

cat("fraction of tiles whose strongest window exceeds a threshold\n")
cat("(max over the four treatments)\n\n")
map_dfr(c(5, 7.5, 10, 15, 20), function(th)
  tibble(threshold = th, n = sum(tl$wald > th),
         `% of tiles` = round(100 * mean(tl$wald > th), 1))) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\ntile size in Mb -- the reason a physical tile is the wrong unit\n\n")
tl %>% summarise(`median Mb` = round(median(Mb), 2),
                    `min Mb` = round(min(Mb), 2),
                    `max Mb` = round(max(Mb), 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nthe straddling tiles, kept:\n\n")
tl %>% filter(f_eu < 1) %>% arrange(chr, tile) %>%
  transmute(chr, cM = paste0(round(cM_lo,1), "-", round(cM_hi,1)),
            Mb = round(Mb, 2), `% euchromatic` = round(100*f_eu),
            wald = round(wald, 1)) %>%
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

# ── peak widths ──────────────────────────────────────────────────────────────
# Well-separated peaks: take the tallest window, exclude +/-5 cM around it,
# repeat. The exclusion is generous on purpose -- the point is to have peaks that
# cannot be shoulders of each other, not to count loci.
cat("\n\nWELL-SEPARATED PEAKS above -log10 P = 15 (euchromatin, +/-5 cM exclusion)\n")
cat("interval = the span within 2 -log10 P of the peak\n\n")
allw <- w %>% filter(is_eu) %>% group_by(chr, pos, cM) %>%
  summarise(wald = max(Wald_log10p), .groups = "drop") %>% arrange(chr, pos)
pk <- list(); x <- allw
while (nrow(x) && max(x$wald) > 15) {
  p <- x[which.max(x$wald), ]; pk[[length(pk) + 1]] <- p
  x <- x %>% filter(!(chr == p$chr & abs(cM - p$cM) <= 5))
}
pk <- bind_rows(pk)
int <- pk %>% pmap_dfr(function(chr, pos, cM, wald) {
  y <- allw %>% filter(chr == !!chr) %>% arrange(pos); k <- which(y$pos == pos)
  hi <- wald - 2
  lo <- k; while (lo > 1  && y$wald[lo - 1] > hi) lo <- lo - 1
  up <- k; while (up < nrow(y) && y$wald[up + 1] > hi) up <- up + 1
  tibble(chr, Mb = round(pos/1e6, 2), wald = round(wald, 1),
         cM_width = round(y$cM[up] - y$cM[lo], 2),
         kb_width = round((y$pos[up] - y$pos[lo]) / 1e3))
})
as.data.frame(int) %>% print(row.names = FALSE)
cat(sprintf("\n  %d peaks;  cM width median %.2f, range %.2f-%.2f\n",
            nrow(int), median(int$cM_width), min(int$cM_width), max(int$cM_width)))
o <- sort(int$cM_width)
cat(sprintf("  %d of %d between %.2f and %.2f cM\n",
            nrow(int) - 1, nrow(int), o[1], o[nrow(int) - 1]))
cat(sprintf("  physical width %.0f kb to %.2f Mb\n",
            min(int$kb_width), max(int$kb_width) / 1e3))

cat("\nchr3L 9.3 Mb by treatment:\n\n")
w %>% filter(chr == "chr3L", abs(pos/1e6 - 9.30) < 0.03) %>%
  group_by(sugar, sex) %>% summarise(wald = round(max(Wald_log10p), 1), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
