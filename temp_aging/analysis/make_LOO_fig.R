# make_LOO_fig.R -- can nine cages predict the tenth?
#
#   Rscript temp_aging/analysis/make_LOO_fig.R
#
# Leave one cage out. From the other nine take e_f = mean(Delta p_f / i) at each
# peak -- the weights. The held-out cage's selected proportion is known, so its
# predicted Delta p is i x e with NOTHING fitted to it. Panel A is those
# predictions against what the cage actually did.
#
# Positions here are the 7 zoom loci, the only per-replicate data held locally.
# The genome-wide version needs the 40 fold-scans on HPC3; the figure is the same.
#
# Panel B is the h2 accounting. A single cage's apparent h2 is inflated by its own
# sampling noise. The cross-product of training and held-out e is unbiased for the
# true h2, since the two carry independent noise, so E[e_tr x e_te] = e_true^2 --
# no s^2/R correction and no winner's curse.

suppressMessages({library(tidyverse); library(patchwork)})
DROP <- c(8,9); BASE <- 13
FOUNDER_COL <- c(AB8="#E69F00", B1="#56B4E9", B2="#009E73", B3="#F0E442",
                 B4="#D55E00", B5="#CC79A7", B6="#000000", B7="#990099")
TRT_LAB <- c(SY10_F="SY10 female", SY20_F="SY20 female",
             SY10_M="SY10 male",   SY20_M="SY20 male")

prop <- map_dfr(c("SY10","SY20"), function(s) map_dfr(c("F","M"), function(x)
  read.table(sprintf("helpfiles/AGE_SY/nov_only/AGE_%s_%s.no89.txt", s, x), header=TRUE) %>%
    filter(TRT=="Z") %>% transmute(sugar=s, sex=x, REP=as.integer(REP), P=Proportion))) %>%
  mutate(i = dnorm(qnorm(1-P))/P)              # selection intensity, known per cage

m <- read.table("process/AGE_SY/AGE_SY_zoom_means.txt.gz", header=TRUE, sep="\t") %>%
  as_tibble() %>% filter(!REP %in% DROP) %>%
  pivot_wider(names_from=TRT, values_from=freq) %>%
  mutate(D = Z - C, chr = sub(":.*","",locus)) %>% filter(pos == peak_pos) %>%
  inner_join(prop, by=c("sugar","sex","REP")) %>%
  mutate(e = D/i, k = ifelse(chr=="chrX" & sex=="M", 1, 2),
         trt = factor(TRT_LAB[paste0(sugar,"_",sex)], levels=TRT_LAB))

# the held-out cage contributes nothing to its own weights
loo <- m %>% group_by(sugar, sex, locus, founder) %>%
  mutate(e_tr = (sum(e)-e)/(n()-1), pbar = mean(C)) %>% ungroup() %>%
  mutate(pred = i * e_tr)

lab <- loo %>% group_by(trt) %>%
  summarise(r = cor(pred, D), slope = coef(lm(D ~ pred))[2], .groups="drop") %>%
  mutate(txt = sprintf("r = %.2f\nslope = %.2f", r, slope))
L <- max(abs(c(loo$pred, loo$D)))

A <- ggplot(loo, aes(pred, D, colour = founder)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey45") +
  geom_point(size=1.5, alpha=0.75) +
  geom_smooth(method="lm", se=FALSE, colour="grey20", linewidth=0.6, formula=y~x) +
  geom_text(data=lab, aes(-L*0.95, L*0.88, label=txt), inherit.aes=FALSE,
            hjust=0, size=3.6, lineheight=0.95) +
  scale_colour_manual(values=FOUNDER_COL, name=NULL) +
  coord_fixed(xlim=c(-L,L), ylim=c(-L,L)) +
  facet_wrap(~trt, nrow=1) +
  labs(x = "predicted change in founder frequency, from the other nine cages",
       y = "observed, held-out cage") +
  guides(colour = guide_legend(nrow=1, override.aes=list(size=3))) +
  theme_bw(BASE) + theme(panel.grid.minor=element_blank(), legend.position="bottom",
                         strip.background=element_rect(fill="grey92", colour=NA))

h2 <- loo %>% group_by(trt, REP) %>%
  summarise(`one cage\nsquared`   = sum(100*k*e^2/pmax(pbar,0.01)),
            `nine cages\nsquared`  = sum(100*k*e_tr^2/pmax(pbar,0.01)),
            `nine x held-out\nunbiased` = sum(100*k*e_tr*e/pmax(pbar,0.01)), .groups="drop") %>%
  pivot_longer(-c(trt,REP), names_to="est", values_to="h2") %>%
  mutate(est = factor(est, levels=c("one cage\nsquared","nine cages\nsquared",
                                    "nine x held-out\nunbiased")))

B <- ggplot(h2, aes(est, h2, colour=est)) +
  geom_hline(yintercept=0, colour="grey80", linewidth=0.3) +
  geom_point(position=position_jitter(width=0.14, height=0), size=1.3, alpha=0.5) +
  stat_summary(fun=mean, geom="crossbar", width=0.55, linewidth=0.4, colour="grey15") +
  scale_colour_manual(values=c("grey60","#1F4E9C","#D62728"), guide="none") +
  facet_wrap(~trt, nrow=1) +
  labs(x=NULL, y=expression(italic(h)^2~"at these positions (%)"),
       subtitle=expression(paste(italic(h)^2,
         " sums squared frequency changes, so noise inflates it. Multiplying two independent estimates instead of squaring one does not."))) +
  theme_bw(BASE) + theme(panel.grid.minor=element_blank(),
                         axis.text.x=element_text(size=BASE-3, lineheight=0.9),
                         plot.subtitle=element_text(size=BASE-2, colour="grey25"),
                         strip.background=element_rect(fill="grey92", colour=NA))

png("temp_aging/analysis/LOO_fig.png", width=12, height=7.2, units="in", res=300)
print(A / B + plot_layout(heights=c(1.05,1)) + plot_annotation(tag_levels="A"))
invisible(dev.off()); cat("Wrote: temp_aging/analysis/LOO_fig.png\n")

cat("\nrho_w (fraction of the position h2 the score captures) and implied individual R2:\n")
print(as.data.frame(h2 %>% group_by(trt, est) %>% summarise(h2=mean(h2), .groups="drop") %>%
  pivot_wider(names_from=est, values_from=h2) %>%
  transmute(trt, rho_w = round(`cross-product`/`mean of nine`, 3),
            R2_individual = round(`cross-product`*`cross-product`/`mean of nine`/100, 4))),
  row.names=FALSE)
