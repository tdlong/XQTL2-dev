# Fit the h2 floor SEPARATELY for each of the four characters, on that
# character's own Wald-null windows, instead of one global F(b) for all four.
suppressMessages(library(tidyverse))
HET<-tribble(~chr,~eu_start,~eu_end,"chrX",2.5,21.2,"chr2L",0.5,22.9,
             "chr2R",1.3,25.1,"chr3L",0.7,24.0,"chr3R",4.5,32.0)
l<-read.table("process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
d<-read.table("process/AGE_SY/AGE_SY_4scan.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
mw<-d%>%group_by(chr,pos)%>%summarise(mwald=max(Wald_log10p),.groups="drop")

# global floor, as published
g<-l%>%group_by(chr,pos)%>%summarise(H2=mean(Cutl_H2),bb=mean(Cutl_H2_bias),
     mwald=max(Wald_log10p),.groups="drop")%>%filter(mwald<2)%>%arrange(bb)
Fg<-approxfun(g$bb,isoreg(g$bb,g$H2)$yf,rule=2)

# per-character floors
mk<-function(dd){dd<-dd%>%arrange(bb); approxfun(dd$bb,isoreg(dd$bb,dd$H2)$yf,rule=2)}
pc<-l%>%inner_join(mw,by=c("chr","pos"))%>%filter(mwald<2)%>%
    group_by(sugar,sex,chr,pos)%>%summarise(H2=mean(Cutl_H2),bb=mean(Cutl_H2_bias),.groups="drop")%>%
    group_by(sugar,sex)%>%group_split()
Fs<-set_names(map(pc,mk), map_chr(pc,~paste0(.x$sugar[1],"_",.x$sex[1])))
cat("floor at the median reported b, per character:\n")
bmed<-median(l$Cutl_H2_bias)
cat(sprintf("  global   %.3f\n",Fg(bmed)))
for(k in names(Fs)) cat(sprintf("  %-9s %.3f\n",k,Fs[[k]](bmed)))

app<-l%>%mutate(key=paste0(sugar,"_",sex),
   v=Cutl_H2-map2_dbl(key,Cutl_H2_bias,~Fs[[.x]](.y)))%>%inner_join(mw,by=c("chr","pos"))
cat("\ncheck: median corrected h2 at Wald-null windows should now be ~0 for all four\n")
app%>%filter(mwald<2)%>%group_by(sugar,sex)%>%
  summarise(med=round(median(v),4),sd=round(sd(v),3),.groups="drop")%>%
  as.data.frame()%>%print(row.names=FALSE)

part<-function(dat){
  w<-dat%>%select(chr,pos,sugar,sex,half,v)%>%
    pivot_wider(names_from=c(sugar,sex,half),values_from=v)%>%drop_na()%>%
    transmute(chr,a=(SY10_F_odd+SY10_F_even)/2,b=(SY20_F_odd+SY20_F_even)/2,
      cc=(SY10_M_odd+SY10_M_even)/2,dd=(SY20_M_odd+SY20_M_even)/2,
      msr=((SY10_F_odd-SY10_F_even)^2+(SY20_F_odd-SY20_F_even)^2+
           (SY10_M_odd-SY10_M_even)^2+(SY20_M_odd-SY20_M_even)^2)/2/4)%>%
    mutate(y=(a+b+cc+dd)/4,F_=(a+b)/2,M_=(cc+dd)/2,S10=(a+cc)/2,S20=(b+dd)/2,
      main=8*y^2-msr,sex=4*((F_-y)^2+(M_-y)^2)-msr,
      diet=4*((S10-y)^2+(S20-y)^2)-msr,
      int=2*((a-F_-S10+y)^2+(b-F_-S20+y)^2+(cc-M_-S10+y)^2+(dd-M_-S20+y)^2)-msr)
  s<-function(z)sum(z,na.rm=TRUE)
  for(sc in c("all","autosomes")){x<-if(sc=="autosomes") w%>%filter(chr!="chrX") else w
    t<-s(x$main)+s(x$sex)+s(x$diet)+s(x$int)
    cat(sprintf("  %-10s main %5.1f  sex %5.1f  diet %4.1f  int %5.1f\n",
      sc,100*s(x$main)/t,100*s(x$sex)/t,100*s(x$diet)/t,100*s(x$int)/t))}
}
cat("\nsplit with GLOBAL floor (as published):\n")
part(l%>%mutate(v=Cutl_H2-Fg(Cutl_H2_bias)))
cat("split with PER-CHARACTER floors:\n"); part(app)

cat("\nper-window M vs F medians, euchromatin, per-character floors:\n")
app%>%group_by(chr,pos,sex)%>%summarise(h2=mean(v),.groups="drop")%>%
  left_join(HET,by="chr")%>%filter(pos/1e6>=eu_start,pos/1e6<=eu_end)%>%
  pivot_wider(names_from=sex,values_from=h2)%>%group_by(chr)%>%
  summarise(`med M`=round(median(M),3),`med F`=round(median(F),3),
    `M>F at`=sprintf("%.0f%%",100*mean(M>F)),.groups="drop")%>%
  as.data.frame()%>%print(row.names=FALSE)
