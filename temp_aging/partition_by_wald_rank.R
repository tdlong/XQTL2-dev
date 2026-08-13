# partition_by_wald_rank.R — how much of the heritability sits in the strongest
# part of the genome, and how much of that is sex + diet.
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript temp_aging/partition_by_wald_rank.R
#
# UNITS. Raw summed window h2 overcounts -- neighbouring windows measure
# overlapping bits of the same signal -- so the genome is rescaled to total 50%,
# the known broad-sense heritability of longevity. That is a real number, not a
# convenience.
#
# BLOCKS are 2 cM and non-overlapping. Heterochromatin is INCLUDED but collapsed:
# each arm contributes exactly one telomeric and one centromeric block, so it
# cannot fragment into dozens of blocks and dominate the ranking. Each of those
# choices matters:
#
#   cM not Mb   h2 is smeared over ~9 Mb in the low-recombination regions
#               flanking the centromeres and ~1 Mb in euchromatin, so a fixed
#               physical block counts one region as ten. A 2 cM block here runs
#               0.38 to 6.38 Mb depending on where it sits.
#   heterochromatin  7% of the sequence but 20% of the h2, a threefold
#               enrichment, and it reads ~0% sex+diet. Split at 2 cM it would
#               contribute many blocks; collapsed it is 8 blocks (4 arms x 2
#               ends), which is what it is worth.
#   no overlap  the scan uses 75 kb windows stepping 5 kb, so windows overlap
#               15-fold. Every 15th window tiles the genome exactly once.
#
# THE INTERVAL integrates two sources of uncertainty in one bootstrap:
#
#   1. which blocks you happened to sample  -- resample 8 cM groups
#   2. the h2 floor                         -- the floor is an isotonic fit of h2
#      on the reported bias over Wald-null windows; resample those windows (in
#      blocks) and refit, so each replicate carries its own floor
#
# The floor is the larger of the two and is the reason low-h2 regions have such
# wide intervals: there main is near zero, so the ratio has a pole in it.

suppressMessages({library(tidyverse); library(splines)})

LONG   <- "process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz"
FLYMAP <- "pipeline/helpfiles/flymap.r6.txt"
CM     <- 2       # block width, cM
WNULL  <- 2       # Wald below this: frequencies did not move, so true h2 ~ 0
NBOOT  <- 300

HET <- tribble(~chr, ~eu_start, ~eu_end,
  "chrX", 2.5, 21.2, "chr2L", 0.5, 22.9, "chr2R", 1.3, 25.1,
  "chr3L", 0.7, 24.0, "chr3R", 4.5, 32.0)

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

# The X is INCLUDED. It contributes heritability and data, and the question this
# table answers is how well % non-main can be estimated -- not whether the X's
# sex component transfers to another species.
long <- read.table(LONG, header = TRUE, sep = "\t") %>% as_tibble()

# windows used to fit the floor: Wald says nothing moved, so the observed h2 IS
# the floor there
nullw <- long %>% group_by(chr, pos) %>%
  summarise(H2 = mean(Cutl_H2), b = mean(Cutl_H2_bias),
            wald = max(Wald_log10p), .groups = "drop") %>%
  filter(wald < WNULL) %>% mutate(nb = paste0(chr, "_", pos %/% 2e6))

# geometry, computed once
geo <- long %>% distinct(chr, pos) %>% add_genetic() %>%
  left_join(HET, by = "chr") %>%
  arrange(chr, pos) %>% group_by(chr) %>% slice(seq(1, n(), by = 15)) %>% ungroup() %>%
  mutate(zone = case_when(pos/1e6 <  eu_start ~ "start",
                          pos/1e6 >  eu_end   ~ "end",
                          TRUE                ~ "eu"),
         # euchromatin: 2 cM blocks. heterochromatin: ONE block per arm end.
         blk = ifelse(zone == "eu", paste0(chr, "_", floor(cM / CM)),
                                    paste0(chr, "_het_", zone)),
         grp = ifelse(zone == "eu", paste0(chr, "_", floor(cM / 8)),
                                    paste0(chr, "_het_", zone)),
         het = zone != "eu")

wide <- long %>% select(chr, pos, sex, sugar, half, Cutl_H2, Cutl_H2_bias, Wald_log10p) %>%
  inner_join(geo %>% select(chr, pos, blk, grp, het), by = c("chr","pos"))

# terms per 2 cM block, given a floor function
blocks_for <- function(bcal) {
  wide %>% mutate(v = Cutl_H2 - bcal(Cutl_H2_bias)) %>%
    select(chr, pos, blk, grp, het, sex, sugar, half, v, Wald_log10p) %>%
    group_by(chr, pos, blk, grp, het) %>%
    mutate(wald = max(Wald_log10p)) %>% ungroup() %>% select(-Wald_log10p) %>%
    pivot_wider(names_from = c(sugar, sex, half), values_from = v) %>% drop_na() %>%
    transmute(blk, grp, het, wald,
      a=(SY10_F_odd+SY10_F_even)/2, b=(SY20_F_odd+SY20_F_even)/2,
      cc=(SY10_M_odd+SY10_M_even)/2, d=(SY20_M_odd+SY20_M_even)/2,
      msr=((SY10_F_odd-SY10_F_even)^2+(SY20_F_odd-SY20_F_even)^2+
           (SY10_M_odd-SY10_M_even)^2+(SY20_M_odd-SY20_M_even)^2)/2/4) %>%
    mutate(H2=(a+b+cc+d)/4, F_=(a+b)/2, M_=(cc+d)/2, S10=(a+cc)/2, S20=(b+d)/2,
      main=8*H2^2-msr, sex=4*((F_-H2)^2+(M_-H2)^2)-msr,
      diet=4*((S10-H2)^2+(S20-H2)^2)-msr,
      int=2*((a-F_-S10+H2)^2+(b-F_-S20+H2)^2+(cc-M_-S10+H2)^2+(d-M_-S20+H2)^2)-msr,
      other=sex+diet+int) %>%
    group_by(blk, grp, het) %>%
    summarise(H2=sum(H2), main=sum(main), other=sum(other), wald=max(wald),
              tiles=n(), Mb=n()*0.075, .groups="drop") %>%
    filter(tiles >= 4 | het) %>%
    mutate(tot=main+other, fo=ifelse(tot>0, other/tot, NA),
           h2_other=H2*fo, h2_main=H2*(1-fo))
}

fit_floor <- function(d) approxfun(d$b, isoreg(d$b, d$H2)$yf, rule = 2)

obs <- blocks_for(fit_floor(nullw %>% arrange(b)))
CUTS <- c(5, 10, 20, 30, 50, 100)
stat <- function(bl) map_dbl(CUTS, function(pc) {
  s <- bl %>% filter(wald >= quantile(bl$wald, max(0, 1 - pc/100)))
  100 * sum(s$h2_other, na.rm = TRUE) / sum(s$H2)
})
share <- map_dbl(CUTS, function(pc) {
  s <- obs %>% filter(wald >= quantile(obs$wald, max(0, 1 - pc/100)))
  100 * sum(s$H2) / sum(obs$H2)
})
nblk <- map_int(CUTS, function(pc)
  nrow(obs %>% filter(wald >= quantile(obs$wald, max(0, 1 - pc/100)))))

set.seed(1)
nb_ids <- unique(nullw$nb); nb_ix <- split(seq_len(nrow(nullw)), nullw$nb)
boot <- matrix(NA_real_, NBOOT, length(CUTS))
for (i in seq_len(NBOOT)) {
  nd <- nullw[unlist(nb_ix[sample(nb_ids, length(nb_ids), TRUE)]), ] %>% arrange(b)
  bl <- blocks_for(fit_floor(nd))
  g <- unique(bl$grp); ix <- split(seq_len(nrow(bl)), bl$grp)
  boot[i, ] <- stat(bl[unlist(ix[sample(g, length(g), TRUE)]), ])
}

H2_TOT <- 50
K <- H2_TOT / sum(obs$H2)
sel <- function(pc) obs %>% filter(wald >= quantile(obs$wald, max(0, 1 - pc/100)))

cat(sprintf("whole genome: %d blocks -- %d euchromatic of %g cM, %d heterochromatic\n",
            nrow(obs), sum(!obs$het), CM, sum(obs$het)))
cat("            (one telomeric + one centromeric block per arm)\n")
cat(sprintf("h2 rescaled so the genome totals %g%%; interval from %d bootstraps over\n",
            H2_TOT, NBOOT))
cat("BOTH the block sample and the floor fit\n\n")

tibble(`top % genome` = CUTS,
       # euchromatic + heterochromatic, so you can see the arm-end blocks enter
       blocks = map_chr(CUTS, ~sprintf("%d+%d", sum(!sel(.x)$het), sum(sel(.x)$het))),
       Mb     = map_dbl(CUTS, ~round(sum(sel(.x)$Mb), 1)),
       `h2 (of 50)` = map_dbl(CUTS, ~round(K * sum(sel(.x)$H2), 1)),
       main   = map_dbl(CUTS, ~round(K * sum(sel(.x)$h2_main, na.rm=TRUE), 1)),
       `sex+diet` = map_dbl(CUTS, ~round(K * sum(sel(.x)$h2_other, na.rm=TRUE), 2)),
       `% sex+diet` = round(stat(obs), 1),
       CI = sprintf("[%.1f, %.1f]", apply(boot, 2, quantile, .025, na.rm = TRUE),
                                    apply(boot, 2, quantile, .975, na.rm = TRUE))) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nthe heterochromatic blocks on their own:\n")
obs %>% filter(het) %>% arrange(desc(H2)) %>%
  transmute(block = blk, Mb = round(Mb, 1), `h2 (of 50)` = round(K * H2, 2),
            wald = round(wald, 1), `% sex+diet` = round(100 * h2_other / H2, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)
