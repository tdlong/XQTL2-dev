suppressMessages(library(tidyverse))
m<-read.table("process/AGE_SY/AGE_SY_zoom_means.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
s<-read.table("process/AGE_SY/AGE_SY_4scan.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
best<-m%>%distinct(locus,chr,peak_pos)%>%pmap_dfr(function(locus,chr,peak_pos)
  s%>%filter(chr==!!chr,pos==!!peak_pos)%>%slice_max(Wald_log10p,n=1)%>%
    transmute(locus,sugar,sex))
d<-m%>%inner_join(best,by=c("locus","sugar","sex"))%>%
  pivot_wider(names_from=TRT,values_from=freq,names_prefix="f")%>%
  filter(pos==peak_pos)
rare<-d%>%group_by(locus,founder)%>%summarise(mC=mean(fC),.groups="drop")%>%filter(mC<0.025)
drop<-d%>%anti_join(rare,by=c("locus","founder"))%>%group_by(locus,founder)%>%
  summarise(dZC=mean(fZ-fC),.groups="drop")%>%group_by(locus)%>%slice_min(dZC,n=1)%>%ungroup()
cat("At each peak: the founder that drops most under selection, and what its\n")
cat("frequency does in the CONTROLS across the 12 replicates (= over time).\n\n")
d%>%inner_join(drop,by=c("locus","founder"))%>%group_by(locus,founder,dZC)%>%
  summarise(ctrl_slope_per_rep=coef(lm(fC~REP))[2],
            p=summary(lm(fC~REP))$coefficients[2,4],
            ctrl_change_over_12=coef(lm(fC~REP))[2]*11,
            sd_ctrl=sd(fC),.groups="drop")%>%
  transmute(locus,founder,
    `drop in selected`=round(dZC,3),
    `control change over 12 reps`=round(ctrl_change_over_12,3),
    `p (slope)`=round(p,3),
    `ratio`=round(ctrl_change_over_12/dZC,2))%>%
  arrange(`drop in selected`)%>%as.data.frame()%>%print(row.names=FALSE)
cat("\nratio near 0 = control flat, drop is selection.\n")
cat("ratio near 1 = the founder is going out of the cage anyway.\n")
