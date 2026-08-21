# Everything restricted to the significant part of the genome: non-overlapping
# euchromatic tiles with Wald > 5 in at least one character.
#
# DATA: the 10-replicate AGE_SY dataset -- replicates 8 and 9 dropped, those
# being the May 2023 cage (helpfiles/AGE_2024/population_assignment.txt). The
# 12-replicate files still exist on disk; nothing here reads them.
suppressMessages(library(tidyverse))


HET<-tribble(~chr,~eu_start,~eu_end,"chrX",2.5,21.2,"chr2L",0.5,22.9,
             "chr2R",1.3,25.1,"chr3L",0.7,24.0,"chr3R",4.5,32.0)
l<-read.table("process/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
d<-read.table("process/AGE_SY/AGE_SY_4scan_no89.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
tiles<-d%>%distinct(chr,pos)%>%arrange(chr,pos)%>%group_by(chr)%>%
       slice(seq(1,n(),by=15))%>%ungroup()
sig<-d%>%group_by(chr,pos)%>%summarise(mw=max(Wald_log10p),.groups="drop")%>%
     inner_join(tiles,by=c("chr","pos"))%>%left_join(HET,by="chr")%>%
     filter(pos/1e6>=eu_start,pos/1e6<=eu_end)
cat(sprintf("euchromatic tiles: %d, of which Wald>5 in some character: %d (%.0f%%)\n\n",
  nrow(sig),sum(sig$mw>5),100*mean(sig$mw>5)))

g<-l%>%group_by(chr,pos)%>%summarise(H2=mean(Cutl_H2),bb=mean(Cutl_H2_bias),
     mw=max(Wald_log10p),.groups="drop")%>%filter(mw<2)%>%arrange(bb)
F_<-approxfun(g$bb,isoreg(g$bb,g$H2)$yf,rule=2)

run<-function(useF,lab){
  w<-l%>%mutate(v=if(useF) Cutl_H2-F_(Cutl_H2_bias) else Cutl_H2)%>%
    inner_join(sig%>%filter(mw>5)%>%select(chr,pos),by=c("chr","pos"))%>%
    select(chr,pos,sugar,sex,half,v)%>%
    pivot_wider(names_from=c(sugar,sex,half),values_from=v)%>%drop_na()%>%
    transmute(chr,pos,a=(SY10_F_odd+SY10_F_even)/2,b=(SY20_F_odd+SY20_F_even)/2,
      cc=(SY10_M_odd+SY10_M_even)/2,dd=(SY20_M_odd+SY20_M_even)/2,
      msr=((SY10_F_odd-SY10_F_even)^2+(SY20_F_odd-SY20_F_even)^2+
           (SY10_M_odd-SY10_M_even)^2+(SY20_M_odd-SY20_M_even)^2)/2/4)%>%
    mutate(y=(a+b+cc+dd)/4,F2=(a+b)/2,M2=(cc+dd)/2,S10=(a+cc)/2,S20=(b+dd)/2,
      main=8*y^2-msr,sex=4*((F2-y)^2+(M2-y)^2)-msr,
      diet=4*((S10-y)^2+(S20-y)^2)-msr,
      int=2*((a-F2-S10+y)^2+(b-F2-S20+y)^2+(cc-M2-S10+y)^2+(dd-M2-S20+y)^2)-msr)
  s<-function(z)sum(z,na.rm=TRUE)
  for(sc in c("all","autosomes")){
    x<-if(sc=="autosomes") w%>%filter(chr!="chrX") else w
    t<-s(x$main)+s(x$sex)+s(x$diet)+s(x$int)
    cat(sprintf("  %-18s %-10s main %5.1f  sex %5.1f  diet %4.1f  int %5.1f\n",
      lab,sc,100*s(x$main)/t,100*s(x$sex)/t,100*s(x$diet)/t,100*s(x$int)/t))}
  invisible(w)
}
cat("partition over SIGNIFICANT tiles only (%):\n")
w<-run(TRUE,"floor removed"); run(FALSE,"floor left in")
cat("\n  ^ if these two agree, the floor is irrelevant here, which is the point\n")

cat("\nmale vs female h2 in significant tiles, by arm (floor removed):\n\n")
l%>%mutate(v=Cutl_H2-F_(Cutl_H2_bias))%>%
  inner_join(sig%>%filter(mw>5)%>%select(chr,pos),by=c("chr","pos"))%>%
  group_by(chr,pos,sex)%>%summarise(h2=mean(v),.groups="drop")%>%
  pivot_wider(names_from=sex,values_from=h2)%>%group_by(chr)%>%
  summarise(tiles=n(),`med M`=round(median(M),2),`med F`=round(median(F),2),
    `M/F`=round(median(M)/median(F),2),`M>F at`=sprintf("%.0f%%",100*mean(M>F)),
    .groups="drop")%>%as.data.frame()%>%print(row.names=FALSE)
cat("\nsame, floor NOT removed (raw), to show it does not matter here:\n\n")
l%>%inner_join(sig%>%filter(mw>5)%>%select(chr,pos),by=c("chr","pos"))%>%
  group_by(chr,pos,sex)%>%summarise(h2=mean(Cutl_H2),.groups="drop")%>%
  pivot_wider(names_from=sex,values_from=h2)%>%group_by(chr)%>%
  summarise(`M>F at`=sprintf("%.0f%%",100*mean(M>F)),.groups="drop")%>%
  as.data.frame()%>%print(row.names=FALSE)
