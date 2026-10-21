# yardstick.R -- marker-based heritability of each trait.  NOT a hold-out quantity.
#
#   Rscript temp_aging/analysis/yardstick.R
#
# Sum of the full ten-replicate heritability estimate over the 1 cM tiles,
# divided by the number of tiles a single locus is counted in.  That number is
# LD: how far a founder haplotype block extends.  It is read off the footprint of
# the one isolated large locus, 3L 9.30, where h2 falls to half by 2 cM and to a
# tenth by 4.5 cM -- integral over peak = 4.4 tiles.  The autocorrelation of h2
# along the genome gives the same.  SE is a jackknife over replicates.
#
# Per-replicate responses come from loo_u.txt.gz only because that file already
# holds them (the held-out column of each fold IS that replicate); nothing about
# the hold-out design enters the estimate.

suppressMessages(library(tidyverse))
TPL <- 4.4
U <- read.table("temp_aging/numbers/loo_u.txt.gz", header=TRUE, sep="\t") %>% as_tibble() %>%
  mutate(trt=paste0(sugar,"_",sex)) %>% select(trt, chr, tile, founder, k, rep=holdout, u=u_te)
full <- function(z) z %>% group_by(trt,chr,tile,founder,k) %>%
  summarise(m=mean(u), v=var(u), n=n(), .groups="drop") %>%
  mutate(h=100*k*(m^2 - v/n)) %>% group_by(trt) %>% summarise(h2=sum(h)/TPL, .groups="drop")
est <- full(U)
jk <- U %>% group_split(trt) %>% map_dfr(function(z){
  th <- map_dbl(unique(z$rep), ~ full(z %>% filter(rep != .x))$h2)
  tibble(trt=z$trt[1], se=sqrt((length(th)-1)/length(th)*sum((th-mean(th))^2))) })
out <- est %>% left_join(jk, by="trt") %>% mutate(tiles_per_locus=TPL)
write.table(out, "temp_aging/numbers/yardstick.txt", row.names=FALSE, quote=FALSE, sep="\t")
print(as.data.frame(out %>% mutate(across(where(is.double), ~round(.x,2)))), row.names=FALSE)
