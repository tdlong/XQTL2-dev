# scan_resolution.R — how much of the genome responds, and how wide the peaks are.
#
#   Rscript temp_aging/scan_resolution.R
#
# DATA: the 10-replicate AGE_SY dataset -- replicates 8 and 9 dropped, those
# being the May 2023 cage (helpfiles/AGE_2024/population_assignment.txt).
#
# Fractions are of the GENETIC map, not the physical one. A physical tile
# over-weights low-recombination regions: they occupy many Mb but few cM, and
# that is exactly where the broad elevated Wald sits, so a Mb-weighted "fraction
# of the genome responding" counts the same evidence many times over. Both are
# printed so the size of that difference is visible.
#
# Peak widths are reported in cM first for the same reason -- mapping resolution
# is genetic, so cM is the width that should be roughly constant across peaks
# while the physical width tracks local recombination rate.

suppressMessages(library(tidyverse))
# Euchromatin boundaries: read from the pipeline rather than hardcoded. These
# define the grey bands in the figures; the analysis scripts also use them to
# restrict to euchromatin.
HET <- read.table("pipeline/helpfiles/het_bounds.txt", header = TRUE,
                  comment.char = "#") %>% as_tibble()
d <- read.table("process/AGE_SY/AGE_SY_4scan_no89.txt.gz",header=TRUE,sep="\t") %>%
  as_tibble() %>% left_join(HET,by="chr") %>%
  filter(pos/1e6>=eu_start, pos/1e6<=eu_end)

# non-overlapping tiles; each carries the cM it spans
t <- d %>% distinct(chr,pos,cM) %>% arrange(chr,pos) %>% group_by(chr) %>%
  slice(seq(1,n(),by=15)) %>% mutate(dcM = c(diff(cM), NA), dbp = c(diff(pos), NA)) %>%
  filter(!is.na(dcM), dcM > 0) %>% ungroup()
w <- d %>% inner_join(t %>% select(chr,pos,dcM,dbp), by=c("chr","pos")) %>%
  group_by(chr,pos,dcM,dbp) %>% summarise(wald=max(Wald_log10p), .groups="drop")

cat("EUCHROMATIC GENOME ABOVE A THRESHOLD -- by cM and, for contrast, by Mb\n")
cat("(max over the four treatments; ", nrow(w), " non-overlapping tiles, ",
    round(sum(w$dcM)), " cM, ", round(sum(w$dbp)/1e6), " Mb)\n\n", sep="")
map_dfr(c(5,10,15), function(th) tibble(
  threshold = th,
  `% of cM`  = round(100*sum(w$dcM[w$wald>th])/sum(w$dcM),1),
  `% of Mb`  = round(100*sum(w$dbp[w$wald>th])/sum(w$dbp),1))) %>%
  as.data.frame() %>% print(row.names=FALSE)

cat("\nPEAKS: top-down, excluding +/-5 cM around each, euchromatin, Wald > 15\n\n")
allw <- d %>% group_by(chr,pos,cM) %>% summarise(wald=max(Wald_log10p),.groups="drop") %>%
  arrange(chr,pos)
pk <- list(); x <- allw
while (nrow(x) && max(x$wald) > 15) {
  p <- x[which.max(x$wald),]; pk[[length(pk)+1]] <- p
  x <- x %>% filter(!(chr==p$chr & abs(cM-p$cM) <= 5))
}
pk <- bind_rows(pk)
int <- pk %>% pmap_dfr(function(chr,pos,cM,wald){
  y <- allw %>% filter(chr==!!chr) %>% arrange(pos); k <- which(y$pos==pos)
  hi <- wald - 2
  lo <- k; while(lo>1 && y$wald[lo-1]>hi) lo <- lo-1
  up <- k; while(up<nrow(y) && y$wald[up+1]>hi) up <- up+1
  tibble(chr, Mb=round(pos/1e6,2), wald=round(wald,1),
         cM_width=round(y$cM[up]-y$cM[lo],2),
         kb_width=round((y$pos[up]-y$pos[lo])/1e3))
})
as.data.frame(int) %>% print(row.names=FALSE)
cat(sprintf("\n  n = %d peaks;  cM width median %.2f, range %.2f-%.2f\n",
    nrow(int), median(int$cM_width), min(int$cM_width), max(int$cM_width)))
cat(sprintf("  physical width range %.0f kb to %.2f Mb\n",
    min(int$kb_width), max(int$kb_width)/1e3))
cat(sprintf("\n  chr3L peak: %.1f  (max over treatments)\n", max(pk$wald)))
d %>% filter(chr=="chr3L", abs(pos/1e6-9.3)<0.1) %>% group_by(sugar,sex) %>%
  summarise(w=round(max(Wald_log10p),1),.groups="drop") %>%
  as.data.frame() %>% print(row.names=FALSE)
