# h2 vs Wald, per character. Within a character the two map essentially 1:1, so
# the value of h2 that is distinguishable from zero can be read straight off the
# Wald-null windows. Raw uncorrected h2, i.e. exactly what Fig 1b plots.
# One point per non-overlapping window (every 15th) per character, euchromatin.
suppressMessages(library(tidyverse))
TRT_LEV<-c("SY10 female","SY20 female","SY10 male","SY20 male")
TRT_COL<-c("SY10 female"="#F49AC2","SY20 female"="#D62728",
           "SY10 male"="#8EC7E8","SY20 male"="#1F4E9C")
HET<-tribble(~chr,~eu_start,~eu_end,"chrX",2.5,21.2,"chr2L",0.5,22.9,
             "chr2R",1.3,25.1,"chr3L",0.7,24.0,"chr3R",4.5,32.0)
d<-read.table("process/AGE_SY/AGE_SY_4scan.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
x<-d%>%inner_join(d%>%distinct(chr,pos)%>%arrange(chr,pos)%>%group_by(chr)%>%
     slice(seq(1,n(),by=15))%>%ungroup(),by=c("chr","pos"))%>%
   left_join(HET,by="chr")%>%filter(pos/1e6>=eu_start,pos/1e6<=eu_end)%>%
   mutate(trt=factor(paste0(sugar," ",ifelse(sex=="F","female","male")),levels=TRT_LEV))

cat("h2 at Wald-null windows (this character's own Wald < 2)\n\n")
nul<-x%>%filter(Wald_log10p<2)%>%group_by(trt)%>%
  summarise(n=n(),median=round(median(Cutl_H2),2),`95th`=round(quantile(Cutl_H2,.95),2),
            `99th`=round(quantile(Cutl_H2,.99),2),.groups="drop")
as.data.frame(nul)%>%print(row.names=FALSE)

cat("\nh2 at a given Wald, per character (median of windows within 10% of it)\n\n")
at<-function(t,w){v<-x%>%filter(trt==t,abs(log(Wald_log10p)-log(w))<.1)%>%pull(Cutl_H2)
  if(length(v)<5) NA else median(v)}
expand_grid(trt=TRT_LEV,wald=c(2,5,10,15,20,30,50))%>%
  mutate(h2=round(map2_dbl(trt,wald,at),2))%>%
  pivot_wider(names_from=trt,values_from=h2)%>%as.data.frame()%>%print(row.names=FALSE)

cat("\nfraction of the euchromatic genome above each character's null 95th pct\n\n")
x%>%left_join(nul%>%select(trt,thr=`95th`),by="trt")%>%group_by(trt)%>%
  summarise(threshold=first(thr),`% of genome above`=round(100*mean(Cutl_H2>thr),1),
            .groups="drop")%>%as.data.frame()%>%print(row.names=FALSE)

png("temp_aging/h2_vs_wald.png",width=5.5,height=3.6,units="in",res=300)
print(ggplot(x,aes(Wald_log10p,Cutl_H2,colour=trt))+
  geom_point(size=.45,alpha=.4)+
  scale_x_continuous(trans="log1p",breaks=c(0,1,3,10,30,100,200))+
  scale_colour_manual(values=TRT_COL,name=NULL)+coord_cartesian(ylim=c(0,5))+
  labs(x=expression(Wald~-log[10]*italic(P)),y=expression(italic(h)^2))+
  theme_bw(9)+theme(legend.position="bottom",panel.grid.minor=element_blank())+
  guides(colour=guide_legend(nrow=1,override.aes=list(size=2,alpha=1))))
dev.off()
cat("\nwrote temp_aging/h2_vs_wald.png\n")
