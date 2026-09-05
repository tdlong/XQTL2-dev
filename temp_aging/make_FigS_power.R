# make_FigS_power.R -- power against experiment size, four designs.
#
#   Rscript temp_aging/make_FigS_power.R
#
# Supplementary figure for the power supplement. Same target throughout: a
# locus explaining 2% of individual phenotypic variance, h2 = 50%.
#
# X axis is the number of individuals PHENOTYPED -- individuals for the outbred
# panel and for X-QTL, lines for the inbred-line and RIL panels. That is the
# quantity that actually costs, and it is what makes the designs comparable.
#
# Solid = up to the size of experiments currently run. Dotted = plausible
# extension. Thresholds are each design's own realistic genome-wide bar:
# 0.05/1e6 for the SNP panels, 0.05/1e4 for the RIL linkage scan, and the
# paper's own LWP 7.5 for X-QTL -- which is far stricter than a Bonferroni over
# 268 tiles would be, so the X-QTL curve is conservative.

suppressMessages(library(tidyverse))

p <- 0.15; q <- 1-p
Vp <- 1; h2 <- 0.5; Vg <- h2*Vp; Ve <- Vp-Vg
QTL <- 0.02; r <- 10
alpha <- sqrt(QTL/(2*p*q))
vloc_line <- p*q*(2*alpha)^2          # 0.04, inbreeding doubles it
vbg_line  <- 2*(Vg-QTL)               # 0.96, and the background too
Vline     <- vloc_line + vbg_line + Ve/r
R2_line   <- vloc_line/Vline          # 0.0381

P_SEL <- 0.05; ibar <- dnorm(qnorm(1-P_SEL))/P_SEL; k <- 2

pow_R2 <- function(n, R2, df, thr) pchisq(thr, df, ncp = n*R2/(1-R2), lower.tail = FALSE)
pow_X  <- function(N, thr) pchisq(thr, 7, ncp = (N*P_SEL)*ibar^2*(100*QTL)/(100*k),
                                  lower.tail = FALSE)

thr1 <- qchisq(0.05/1e6, 1, lower.tail = FALSE)
thr7 <- qchisq(0.05/1e4, 7, lower.tail = FALSE)
thrX <- qchisq(10^-7.5,  7, lower.tail = FALSE)

grid <- function(a,b) exp(seq(log(a), log(b), length.out = 300))
mk <- function(lab, n, pw, solid_to) tibble(design = lab, n = n, power = pw,
                                            part = ifelse(n <= solid_to, "run", "extension"))

d <- bind_rows(
  mk("Outbred panel, 1 measure",   grid(100, 3000),   pow_R2(grid(100,3000),  QTL/Vp,  1, thr1), 1000),
  mk("Inbred lines, 10 measures",  grid(50,  2000),   pow_R2(grid(50,2000),   R2_line, 1, thr1), 200),
  mk("8-way RILs, 10 measures",    grid(100, 1500),   pow_R2(grid(100,1500),  R2_line, 7, thr7), 750),
  mk("X-QTL (this experiment)",    grid(1000, 2e5),   pow_X(grid(1000,2e5),   thrX),             114000))

pts <- tribble(~design, ~n, ~power,
  "Outbred panel, 1 measure",  1000,  pow_R2(1000, QTL/Vp,  1, thr1),
  "Inbred lines, 10 measures",  200,  pow_R2(200,  R2_line, 1, thr1),
  "8-way RILs, 10 measures",    750,  pow_R2(750,  R2_line, 7, thr7),
  "X-QTL (this experiment)", 114000,  pow_X(114000, thrX))
cat("\npower at the size each design is actually run at:\n"); print(as.data.frame(pts), row.names = FALSE)

LEV <- c("X-QTL (this experiment)","8-way RILs, 10 measures",
         "Outbred panel, 1 measure","Inbred lines, 10 measures")
d   <- d   %>% mutate(design = factor(design, levels = LEV))
pts <- pts %>% mutate(design = factor(design, levels = LEV))
COL <- c("#B4413F","#3D6E9C","#4E8A5B","#8A6BA8"); names(COL) <- LEV

g <- ggplot(d, aes(n, power, colour = design)) +
  geom_hline(yintercept = 0.8, linewidth = 0.3, colour = "grey65", linetype = "22") +
  geom_line(aes(linetype = part), linewidth = 0.8) +
  geom_point(data = pts, size = 2.4) +
  scale_x_log10(breaks = c(100,300,1000,3000,10000,30000,100000),
                labels = c("100","300","1,000","3,000","10,000","30,000","100,000")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = expansion(0.01)) +
  scale_colour_manual(values = COL, name = NULL) +
  scale_linetype_manual(values = c(run = "solid", extension = "22"), guide = "none") +
  labs(x = "individuals phenotyped (lines, for the inbred and RIL panels)",
       y = "power to detect a locus explaining 2% of phenotypic variance") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.02, 0.98), legend.justification = c(0,1),
        legend.background = element_rect(fill = alpha("white",0.85), colour = NA),
        legend.key.height = unit(11,"pt"))

ggsave("figures/FigureS_power.png", g, width = 6.6, height = 4.4, dpi = 300)
cat("Wrote figures/FigureS_power.png\n")
cat(sprintf("\nthresholds: SNP panels %.1f (1 df), RILs %.1f (7 df), X-QTL %.1f (7 df)\n",
            thr1, thr7, thrX))
