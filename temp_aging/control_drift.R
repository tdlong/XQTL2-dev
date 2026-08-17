# Are the control trajectories doing anything beyond drift?
suppressMessages(library(tidyverse)); set.seed(1)
m<-read.table("process/AGE_SY/AGE_SY_zoom_means.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
s<-read.table("process/AGE_SY/AGE_SY_4scan.txt.gz",header=TRUE,sep="\t")%>%as_tibble()
L<-m%>%distinct(locus,chr,peak_pos)
b<-L%>%pmap_dfr(function(locus,chr,peak_pos) s%>%filter(chr==!!chr,pos==!!peak_pos)%>%
  slice_max(Wald_log10p,n=1)%>%transmute(locus,sugar,sex))
w<-m%>%inner_join(b,by=c("locus","sugar","sex"))%>%filter(pos==peak_pos)%>%
  pivot_wider(names_from=TRT,values_from=freq,names_prefix="f")%>%
  group_by(locus,founder)%>%mutate(rare=mean(fC)<0.025)%>%ungroup()%>%filter(!rare)
mv<-w%>%group_by(locus,founder)%>%summarise(D=mean(fZ-fC),.groups="drop")%>%
  group_by(locus)%>%{bind_rows(slice_max(.,D,n=1)%>%mutate(dir="protective"),
                               slice_min(.,D,n=1)%>%mutate(dir="susceptible"))}%>%ungroup()
d<-w%>%inner_join(mv%>%select(locus,founder,dir),by=c("locus","founder"))

fit<-d%>%group_by(dir,locus,founder)%>%
  summarise(slope=coef(lm(fC~REP))[2], p=summary(lm(fC~REP))$coefficients[2,4],
            .groups="drop")
cat("slope of control frequency on replicate, 14 series\n\n")
fit%>%mutate(slope=round(slope,4),p=round(p,3))%>%arrange(dir,p)%>%
  as.data.frame()%>%print(row.names=FALSE)
cat(sprintf("\nnominally p<0.05: %d of %d (0.7 expected if every series were flat noise)\n",
    sum(fit$p<0.05), nrow(fit)))

# But these are frequency trajectories, not independent points. A pure random
# walk produces apparent linear trends far more often than 5% of the time.
# Calibrate: simulate drift with the same per-step variance as the observed
# series, and see how often a slope this large appears.
step_sd<-d%>%arrange(locus,founder,REP)%>%group_by(locus,founder)%>%
  summarise(sd=sd(diff(fC)),.groups="drop")
sim_p<-map2_dbl(fit$slope, step_sd$sd, function(obs, sdev){
  z<-replicate(4000, {x<-cumsum(c(0,rnorm(11,0,sdev))); abs(coef(lm(x~seq_along(x)))[2])})
  mean(z >= abs(obs))})
cat("\nsame slopes against a random walk of the same step size:\n\n")
fit%>%mutate(slope=round(slope,4), `p vs iid`=round(p,3),
             `p vs drift`=round(sim_p,3))%>%arrange(dir,`p vs drift`)%>%
  as.data.frame()%>%print(row.names=FALSE)
cat(sprintf("\np<0.05 against drift: %d of %d\n", sum(sim_p<0.05), length(sim_p)))
cat(sprintf("mean slope, protective %.4f  susceptible %.4f  (difference %.4f)\n",
  mean(fit$slope[fit$dir=="protective"]), mean(fit$slope[fit$dir=="susceptible"]),
  mean(fit$slope[fit$dir=="protective"])-mean(fit$slope[fit$dir=="susceptible"])))
cat(sprintf("paired t on the 7 loci, protective vs susceptible slope: p = %.3f\n",
  t.test(fit$slope[fit$dir=="protective"], fit$slope[fit$dir=="susceptible"],
         paired=TRUE)$p.value))
