suppressMessages(library(tidyverse))
FLYMAP<-"pipeline/helpfiles/flymap.r6.txt"
add_genetic<-function(df){fm<-read.table(FLYMAP,header=FALSE);colnames(fm)<-c("chr","pos","cM")
 df$cM<-NA_real_;for(c_ in unique(df$chr)){x<-fm%>%filter(chr==c_)
 o<-ksmooth(x$pos,x$cM,kernel="normal",bandwidth=3e6);df$cM[df$chr==c_]<-splinefun(o$x,o$y)(df$pos[df$chr==c_])};df}
HET<-tribble(~chr,~eu_start,~eu_end,"chrX",2.5,21.2,"chr2L",0.5,22.9,
             "chr2R",1.3,25.1,"chr3L",0.7,24.0,"chr3R",4.5,32.0)
d<-read.table("process/AGE_SY/AGE_SY_4scan.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
w<-d%>%group_by(chr,pos)%>%summarise(wald=max(Wald_log10p),.groups="drop")%>%
   add_genetic()%>%left_join(HET,by="chr")%>%
   mutate(eu=pos/1e6>=eu_start&pos/1e6<=eu_end)%>%arrange(chr,pos)
cat("fraction of the EUCHROMATIC genome above a threshold (max over 4 scans)\n")
for(TH in c(5,10,15,20)) cat(sprintf("  -log10P > %2d : %4.1f%%\n",TH,100*mean(w$wald[w$eu]>TH)))
cat("\nsame, within a single scan\n")
d%>%left_join(HET,by="chr")%>%filter(pos/1e6>=eu_start,pos/1e6<=eu_end)%>%group_by(sugar,sex)%>%
  summarise(`>5`=round(100*mean(Wald_log10p>5),1),`>10`=round(100*mean(Wald_log10p>10),1),
            `>15`=round(100*mean(Wald_log10p>15),1),.groups="drop")%>%
  as.data.frame()%>%print(row.names=FALSE)

# the better peaks: top-down, but exclude generously (5 cM) so shoulders don't count
pk<-list(); x<-w%>%filter(eu)
while(nrow(x) && max(x$wald)>15){p<-x[which.max(x$wald),];pk[[length(pk)+1]]<-p
  x<-x%>%filter(!(chr==p$chr & abs(cM-p$cM)<=5))}
pk<-bind_rows(pk)
cat(sprintf("\ndistinct euchromatic peaks above 15, +/-5 cM exclusion: %d\n\n",nrow(pk)))
int<-function(chr_,pos_,drop){x<-w%>%filter(chr==chr_)%>%arrange(pos);k<-which(x$pos==pos_)
  hi<-x$wald[k]-drop;lo<-k;while(lo>1&&x$wald[lo-1]>hi)lo<-lo-1
  up<-k;while(up<nrow(x)&&x$wald[up+1]>hi)up<-up+1
  c(kb=(x$pos[up]-x$pos[lo])/1e3, cM=x$cM[up]-x$cM[lo])}
pk%>%pmap_dfr(function(chr,pos,wald,cM,...){
  a<-int(chr,pos,2); b<-int(chr,pos,wald/2)
  tibble(chr,Mb=round(pos/1e6,2),wald=round(wald,1),
         `2-unit drop kb`=round(a["kb"]),`2-unit cM`=round(a["cM"],2),
         `half-max kb`=round(b["kb"]),`half-max cM`=round(b["cM"],2))})%>%
  as.data.frame()%>%print(row.names=FALSE)
