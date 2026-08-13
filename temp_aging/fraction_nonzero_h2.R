suppressMessages(library(tidyverse))
set.seed(1)
# Under H0 the main term is 1 df and the replicate error 4 df, both scaled by the
# same sigma^2, so main - MS_rep > 0 happens well below half the time.
p0 <- mean(rchisq(2e6,1) > rchisq(2e6,4)/4)
cat(sprintf("null P(main exceeds replicate error) = %.3f  (NOT 0.5)\n\n", p0))

HET<-tribble(~chr,~eu_start,~eu_end,"chrX",2.5,21.2,"chr2L",0.5,22.9,
             "chr2R",1.3,25.1,"chr3L",0.7,24.0,"chr3R",4.5,32.0)
l<-read.table("process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
nw<-l%>%group_by(chr,pos)%>%summarise(H2=mean(Cutl_H2),b=mean(Cutl_H2_bias),
      wald=max(Wald_log10p),.groups="drop")%>%filter(wald<2)%>%arrange(b)
f<-approxfun(nw$b,isoreg(nw$b,nw$H2)$yf,rule=2)

# every 15th window: tile the genome once, no overlap
tiles <- l%>%distinct(chr,pos)%>%arrange(chr,pos)%>%group_by(chr)%>%
  slice(seq(1,n(),by=15))%>%ungroup()

d<-l%>%inner_join(tiles,by=c("chr","pos"))%>%mutate(v=Cutl_H2-f(Cutl_H2_bias))%>%
  select(chr,pos,sex,sugar,half,v)%>%
  pivot_wider(names_from=c(sugar,sex,half),values_from=v)%>%drop_na()%>%
  transmute(chr,pos,a=(SY10_F_odd+SY10_F_even)/2,b=(SY20_F_odd+SY20_F_even)/2,
    cc=(SY10_M_odd+SY10_M_even)/2,d=(SY20_M_odd+SY20_M_even)/2,
    msr=((SY10_F_odd-SY10_F_even)^2+(SY20_F_odd-SY20_F_even)^2+
         (SY10_M_odd-SY10_M_even)^2+(SY20_M_odd-SY20_M_even)^2)/2/4)%>%
  mutate(H2=(a+b+cc+d)/4, main=8*H2^2-msr)%>%
  left_join(HET,by="chr")%>%mutate(eu=pos/1e6>=eu_start&pos/1e6<=eu_end)

for (lab in c("euchromatin","whole genome")) {
  x <- if (lab=="euchromatin") d%>%filter(eu) else d
  obs <- mean(x$main>0)
  pi0 <- (1-obs)/(1-p0)
  cat(sprintf("%-13s  %d non-overlapping tiles\n", lab, nrow(x)))
  cat(sprintf("   main exceeds replicate error at %.1f%% of tiles (null expects %.1f%%)\n",
      100*obs, 100*p0))
  cat(sprintf("   => at most %.0f%% of tiles are null, i.e. at least %.0f%% carry h2 > 0\n\n",
      100*pi0, 100*(1-pi0)))
}
cat("same, per treatment group (euchromatin, using each group's own two halves):\n")
l%>%inner_join(tiles,by=c("chr","pos"))%>%mutate(v=Cutl_H2-f(Cutl_H2_bias))%>%
  select(chr,pos,sex,sugar,half,v)%>%pivot_wider(names_from=half,values_from=v)%>%
  drop_na()%>%left_join(HET,by="chr")%>%filter(pos/1e6>=eu_start,pos/1e6<=eu_end)%>%
  mutate(m=((odd+even)/2)^2*2 - (odd-even)^2/2)%>%   # 1 df each, mean vs within
  group_by(sugar,sex)%>%summarise(`% > 0`=round(100*mean(m>0),1),.groups="drop")%>%
  as.data.frame()%>%print(row.names=FALSE)
