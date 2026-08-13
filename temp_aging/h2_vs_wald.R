# h2 vs Wald, to read off what value of h2 is distinguishable from zero.
# Two panels, matching what Figure 1 actually plots:
#   B  raw Cutl_H2 from the full 12-replicate scans, UNCORRECTED, as in Fig 1b
#   C  the main-effect component, floor-corrected, in h2 units, as in Fig 1c
# One point per non-overlapping window (every 15th) per character.
suppressMessages(library(tidyverse))
TRT_LEV<-c("SY10 female","SY20 female","SY10 male","SY20 male")
TRT_COL<-c("SY10 female"="#F49AC2","SY20 female"="#D62728",
           "SY10 male"="#8EC7E8","SY20 male"="#1F4E9C")
HET<-tribble(~chr,~eu_start,~eu_end,"chrX",2.5,21.2,"chr2L",0.5,22.9,
             "chr2R",1.3,25.1,"chr3L",0.7,24.0,"chr3R",4.5,32.0)
eu<-function(d) d%>%left_join(HET,by="chr")%>%filter(pos/1e6>=eu_start,pos/1e6<=eu_end)
B<-c(-Inf,2,3,5,7,10,15,20,30,Inf); L<-c("<2","2-3","3-5","5-7","7-10","10-15","15-20","20-30",">30")

## ---- panel B: raw h2 from the full scans -------------------------------
d<-read.table("process/AGE_SY/AGE_SY_4scan.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
tiles<-d%>%distinct(chr,pos)%>%arrange(chr,pos)%>%group_by(chr)%>%
       slice(seq(1,n(),by=15))%>%ungroup()
b<-d%>%inner_join(tiles,by=c("chr","pos"))%>%eu()%>%
   mutate(trt=factor(paste0(sugar," ",ifelse(sex=="F","female","male")),levels=TRT_LEV),
          h2=Cutl_H2)
nb<-b%>%filter(Wald_log10p<2)
TB<-quantile(nb$h2,.95)
cat("=== PANEL B: raw (uncorrected) h2 ===\n")
cat(sprintf("null windows, this character's Wald<2: n=%d\n",nrow(nb)))
cat(sprintf("  raw h2 there: median %.2f, 95th pct %.2f, 99th pct %.2f\n",
    median(nb$h2),TB,quantile(nb$h2,.99)))
cat("  per character (median / 95th):\n")
nb%>%group_by(trt)%>%summarise(n=n(),med=round(median(h2),2),
  `95th`=round(quantile(h2,.95),2),.groups="drop")%>%as.data.frame()%>%print(row.names=FALSE)
cat("\nmedian raw h2 by Wald bin:\n")
b%>%mutate(bin=cut(Wald_log10p,B,labels=L))%>%group_by(bin)%>%
  summarise(n=n(),`median h2`=round(median(h2),2),
    `% > null 95th`=round(100*mean(h2>TB),1),.groups="drop")%>%
  as.data.frame()%>%print(row.names=FALSE)
cat(sprintf("\n  => raw h2 below ~%.2f is not distinguishable from zero.\n",TB))
cat(sprintf("     %.1f%% of euchromatic tiles x characters exceed it.\n",100*mean(b$h2>TB)))
b%>%group_by(trt)%>%summarise(`% above`=round(100*mean(h2>TB),1),.groups="drop")%>%
  as.data.frame()%>%print(row.names=FALSE)

## ---- panel C: main effect, corrected ------------------------------------
l<-read.table("process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
nw<-l%>%group_by(chr,pos)%>%summarise(H2=mean(Cutl_H2),bb=mean(Cutl_H2_bias),
      mw=max(Wald_log10p),.groups="drop")%>%filter(mw<2)%>%arrange(bb)
f<-approxfun(nw$bb,isoreg(nw$bb,nw$H2)$yf,rule=2)
mw<-d%>%group_by(chr,pos)%>%summarise(wald=max(Wald_log10p),.groups="drop")
cc<-l%>%mutate(v=Cutl_H2-f(Cutl_H2_bias))%>%select(chr,pos,sugar,sex,half,v)%>%
  pivot_wider(names_from=c(sugar,sex,half),values_from=v)%>%drop_na()%>%
  transmute(chr,pos,a=(SY10_F_odd+SY10_F_even)/2,b2=(SY20_F_odd+SY20_F_even)/2,
    c2=(SY10_M_odd+SY10_M_even)/2,d2=(SY20_M_odd+SY20_M_even)/2,
    msr=((SY10_F_odd-SY10_F_even)^2+(SY20_F_odd-SY20_F_even)^2+
         (SY10_M_odd-SY10_M_even)^2+(SY20_M_odd-SY20_M_even)^2)/2/4)%>%
  mutate(ybar=(a+b2+c2+d2)/4, ss=8*ybar^2-msr,
         main=sign(ybar)*sqrt(abs(ss)/8))%>%     # h2 units, as plotted in Fig 1c
  inner_join(mw,by=c("chr","pos"))%>%inner_join(tiles,by=c("chr","pos"))%>%eu()
nc<-cc%>%filter(wald<2); TC<-quantile(nc$main,.95)
cat("\n=== PANEL C: main effect, floor-corrected, h2 units ===\n")
cat(sprintf("null windows, max Wald<2: n=%d; main there: median %.2f, 95th %.2f, 99th %.2f\n",
    nrow(nc),median(nc$main),TC,quantile(nc$main,.99)))
cat("\nmedian main by Wald bin:\n")
cc%>%mutate(bin=cut(wald,B,labels=L))%>%group_by(bin)%>%
  summarise(n=n(),`median main`=round(median(main),2),
    `% > null 95th`=round(100*mean(main>TC),1),.groups="drop")%>%
  as.data.frame()%>%print(row.names=FALSE)
cat(sprintf("\n  => main below ~%.2f is not distinguishable from zero.\n",TC))
cat(sprintf("     %.1f%% of euchromatic tiles exceed it.\n",100*mean(cc$main>TC)))

## ---- figure --------------------------------------------------------------
suppressMessages(library(patchwork))
lx<-scale_x_continuous(trans="log1p",breaks=c(0,1,3,10,30,100,200))
p1<-ggplot(b,aes(Wald_log10p,h2,colour=trt))+
  geom_hline(yintercept=TB,linetype="dashed",colour="grey30",linewidth=.3)+
  geom_point(size=.45,alpha=.4)+lx+scale_colour_manual(values=TRT_COL,name=NULL)+
  coord_cartesian(ylim=c(0,5))+labs(x=NULL,y=expression(italic(h)^2~(Fig.~1*b)))+
  theme_bw(9)+theme(legend.position="none",panel.grid.minor=element_blank())
p2<-ggplot(cc,aes(wald,main))+
  geom_hline(yintercept=c(0,TC),linetype=c("solid","dashed"),
             colour=c("grey70","grey30"),linewidth=.3)+
  geom_point(size=.45,alpha=.4,colour="grey20")+lx+
  coord_cartesian(ylim=c(-1,5))+
  labs(x=expression(Wald~-log[10]*italic(P)),y=expression(main~(Fig.~1*c)))+
  theme_bw(9)+theme(panel.grid.minor=element_blank())
lg<-cowplot::get_legend(p1+theme(legend.position="bottom")+
     guides(colour=guide_legend(nrow=1,override.aes=list(size=2,alpha=1))))
png("temp_aging/h2_vs_wald.png",width=5.5,height=5.2,units="in",res=300)
print(p1/p2/lg+plot_layout(heights=c(1,1,.12)))
dev.off()
cat("\nwrote temp_aging/h2_vs_wald.png\n")
