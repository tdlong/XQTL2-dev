suppressMessages(library(tidyverse))
FLYMAP <- "pipeline/helpfiles/flymap.r6.txt"
add_genetic <- function(df){ fm<-read.table(FLYMAP,header=FALSE); colnames(fm)<-c("chr","pos","cM")
  df$cM<-NA_real_; for(c_ in unique(df$chr)){x<-fm%>%filter(chr==c_)
  o<-ksmooth(x$pos,x$cM,kernel="normal",bandwidth=3e6); df$cM[df$chr==c_]<-splinefun(o$x,o$y)(df$pos[df$chr==c_])}; df}
HET <- tribble(~chr,~eu_start,~eu_end,"chrX",2.5,21.2,"chr2L",0.5,22.9,
               "chr2R",1.3,25.1,"chr3L",0.7,24.0,"chr3R",4.5,32.0)

d <- read.table("process/AGE_SY/AGE_SY_4scan.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
w <- d %>% group_by(chr,pos) %>% summarise(wald=max(Wald_log10p), h2=max(Cutl_H2),.groups="drop") %>%
     add_genetic() %>% left_join(HET,by="chr") %>%
     mutate(eu = pos/1e6>=eu_start & pos/1e6<=eu_end) %>% arrange(chr,pos)

TH <- 6
runs <- w %>% filter(eu) %>% group_by(chr) %>%
  mutate(above = wald>TH, grp = cumsum(above != lag(above, default=FALSE))) %>%
  filter(above) %>% group_by(chr,grp) %>%
  summarise(kb=(max(pos)-min(pos))/1e3, cM=max(cM)-min(cM), peak=max(wald),
            at=pos[which.max(wald)]/1e6, .groups="drop")
cat(sprintf("contiguous euchromatic runs above -log10P=%g: %d\n", TH, nrow(runs)))
cat(sprintf("run width kb: median %.0f, quartiles %.0f-%.0f, max %.0f\n",
    median(runs$kb), quantile(runs$kb,.25), quantile(runs$kb,.75), max(runs$kb)))
cat(sprintf("run width cM: median %.2f, quartiles %.2f-%.2f, max %.2f\n\n",
    median(runs$cM), quantile(runs$cM,.25), quantile(runs$cM,.75), max(runs$cM)))
cat("run peak height distribution:\n")
print(round(quantile(runs$peak, c(0,.25,.5,.75,.9,1)),1))

cat("\nlargest 12 runs by peak height (kb = full run, not half-max):\n")
runs %>% slice_max(peak,n=12) %>% transmute(chr, peak_Mb=round(at,2),
  wald=round(peak,1), run_kb=round(kb), run_cM=round(cM,2)) %>%
  as.data.frame() %>% print(row.names=FALSE)

cat("\nmax raw Cutl_H2 at the top runs (uncorrected, % of phenotypic variance):\n")
runs %>% slice_max(peak,n=12) %>% pmap_dfr(function(chr,grp,kb,cM,peak,at)
  w %>% filter(chr==!!chr, abs(pos/1e6-at)<0.05) %>% slice_max(wald,n=1) %>%
    transmute(chr, Mb=round(pos/1e6,2), wald=round(wald,1), h2=round(h2,2))) %>%
  as.data.frame() %>% print(row.names=FALSE)

# Also computed for the draft, from the same inputs:
#   - euchromatic fraction above a range of thresholds, whole and per treatment
#   - selection counts from helpfiles/AGE_SY/info.AGE_SY.txt (40 of 48 cages have
#     them; one row is space- not tab-separated, hence the regex parse)
#   - corrected peak h2, which needs the splithalf floor and so is computed in
#     varcomp_H2.R / partition_by_wald_rank.R rather than here
