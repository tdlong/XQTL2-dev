# What value of h2 in Fig 1b is significant?
#
# Within a character, Wald and h2 trace a single curve -- both are driven by the
# same founder frequency shifts -- so a significance threshold on the Wald scale
# transfers directly onto the h2 scale. We fit h2 on log(Wald) per character and
# read off the fit at -log10 P = 5, the cutoff used throughout.
#
# Raw uncorrected h2, i.e. exactly what Fig 1b plots. One point per
# non-overlapping window (every 15th) per character, euchromatin only.
#
#   Rscript temp_aging/h2_threshold.R
suppressMessages({library(tidyverse); library(mgcv)})
HET<-tribble(~chr,~eu_start,~eu_end,"chrX",2.5,21.2,"chr2L",0.5,22.9,
             "chr2R",1.3,25.1,"chr3L",0.7,24.0,"chr3R",4.5,32.0)
d<-read.table("process/AGE_SY/AGE_SY_4scan.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
x<-d%>%inner_join(d%>%distinct(chr,pos)%>%arrange(chr,pos)%>%group_by(chr)%>%
     slice(seq(1,n(),by=15))%>%ungroup(),by=c("chr","pos"))%>%
   left_join(HET,by="chr")%>%filter(pos/1e6>=eu_start,pos/1e6<=eu_end)%>%
   mutate(trt=paste0(sugar,"_",sex), lw=log(Wald_log10p))
map_dfr(c("SY10_F","SY20_F","SY10_M","SY20_M"), function(t){
  z<-x%>%filter(trt==t); m<-gam(Cutl_H2~s(lw), data=z)
  p<-predict(m,newdata=tibble(lw=log(c(2,5,10))),se.fit=TRUE)
  h5<-as.numeric(p$fit[2])
  tibble(character=t, windows=nrow(z),
    `h2 at Wald 2`=round(as.numeric(p$fit[1]),2),
    `h2 at Wald 5`=round(h5,2), se=round(as.numeric(p$se.fit[2]),3),
    `h2 at Wald 10`=round(as.numeric(p$fit[3]),2),
    `% euchr above h2(Wald 5)`=round(100*mean(z$Cutl_H2>h5),1),
    `% euchr Wald>5`=round(100*mean(z$Wald_log10p>5),1))
})%>%as.data.frame()%>%print(row.names=FALSE)
cat("\nnull windows (this character's own Wald < 2), for the floor:\n")
x%>%filter(Wald_log10p<2)%>%group_by(trt)%>%
  summarise(n=n(),`median h2`=round(median(Cutl_H2),2),
            `95th`=round(quantile(Cutl_H2,.95),2),.groups="drop")%>%
  as.data.frame()%>%print(row.names=FALSE)
