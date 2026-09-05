# make_FigS2.R -- power against experiment size, four designs.
#
#   Rscript temp_aging/make_FigS2.R
#
# THE POINT. X-QTL is much LESS efficient per fly than a design that genotypes
# and phenotypes individuals -- pooling and sequencing only the selected tail
# discards what a replicated line mean keeps. Its advantage is that screening
# flies for longevity is cheap and effectively unlimited, so the experiment can
# be two orders of magnitude larger. Each curve therefore stops where that
# design realistically stops:
#
#   solid  = the size these experiments are actually run at
#   dotted = plausible extension
#     inbred lines  200 -> 2,000   rumoured extension of the panel
#     outbred      1000 -> 3,000   a guess; nothing supports a specific ceiling
#     8-way RILs   no extension    the panel is not expected to exceed 750
#     X-QTL        no extension    114,000 is this experiment, per treatment
#
# Two panels: a locus explaining 2% of individual phenotypic variance, and one
# explaining 0.5% -- the size this experiment actually resolves.
#
# Base population: Vp = 1, h2 = 50%, MAF 0.15 for the two SNP designs. Inbred
# lines and RILs are
# phenotyped 10 times each; inbreeding doubles both the locus variance and the
# unlinked background so those cancel, and the gain is replication cutting Ve.
#
# Thresholds are each design's own genome-wide bar: 0.05/1e6 for the SNP panels,
# 0.05/1e4 for the RIL linkage scan, and the paper's own LWP 7.5 for X-QTL --
# stricter than a Bonferroni over 268 tiles, so X-QTL is handicapped here.

suppressMessages(library(tidyverse))

OUT  <- "temp_aging/FigureS2_plot.png"
LOCI <- c("locus explains 2% of phenotypic variance"    = 0.02,
          "locus explains 0.5%, the size resolved here" = 0.005)
p    <- 0.15; q <- 1-p
Vp   <- 1; h2 <- 0.5; Vg <- h2*Vp; Ve <- Vp-Vg; r <- 10
P_SEL <- 0.05; ibar <- dnorm(qnorm(1-P_SEL))/P_SEL; k <- 2

thr1 <- qchisq(0.05/1e6, 1, lower.tail = FALSE)
thr7 <- qchisq(0.05/1e4, 7, lower.tail = FALSE)
thrX <- qchisq(10^-7.5,  7, lower.tail = FALSE)

# inbreeding doubles the locus variance AND the unlinked background, so those
# cancel in the ratio; the gain is replication cutting Ve
R2_line <- function(Q) { vl <- 2*Q; vb <- 2*(Vg-Q); vl/(vl + vb + Ve/r) }

pow_R2 <- function(n, R2, df, thr) pchisq(thr, df, ncp = n*R2/(1-R2), lower.tail = FALSE)
pow_X  <- function(N, Q, thr) pchisq(thr, 7, ncp = (N*P_SEL)*ibar^2*(100*Q)/(100*k),
                                     lower.tail = FALSE)
curve_for <- function(design, n, Q) switch(design,
  "Outbred panel, 1 measure"  = pow_R2(n, Q/Vp,       1, thr1),
  "Inbred lines, 10 measures" = pow_R2(n, R2_line(Q), 1, thr1),
  "8-way RILs, 10 measures"   = pow_R2(n, R2_line(Q), 7, thr7),
  "X-QTL (this experiment)"   = pow_X(n, Q, thrX))

SPEC <- tribble(
  ~design,                       ~lo,  ~run,   ~reach,
  "Outbred panel, 1 measure",     100,   1000,   3000,
  "Inbred lines, 10 measures",     50,    200,   2000,
  "8-way RILs, 10 measures",      100,    750,    750,
  "X-QTL (this experiment)",     2000, 114000, 114000)

grid <- function(a, b) exp(seq(log(a), log(b), length.out = 400))
d <- map_dfr(names(LOCI), function(lab) {
  Q <- LOCI[[lab]]
  pmap_dfr(SPEC, function(design, lo, run, reach) {
    nn <- sort(unique(c(grid(lo, reach), run)))
    pw <- curve_for(design, nn, Q)              # compute before tibble(), or the
    tibble(locus = lab, design = design, n = nn, power = pw,  # column masks the arg
           part = ifelse(nn <= run, "run", "reach"))
  })
})
pts <- map_dfr(names(LOCI), function(lab) {
  Q <- LOCI[[lab]]
  SPEC %>% transmute(locus = lab, design, n = run,
                     power = map2_dbl(design, run, ~curve_for(.x, .y, Q)))
})

cat("power at the size each design is actually run at:\n")
print(as.data.frame(pts %>% mutate(power = round(power, 4))), row.names = FALSE)

LEV <- c("X-QTL (this experiment)", "8-way RILs, 10 measures",
         "Outbred panel, 1 measure", "Inbred lines, 10 measures")
COL <- c("#B4413F", "#3D6E9C", "#4E8A5B", "#8A6BA8"); names(COL) <- LEV
d   <- d   %>% mutate(design = factor(design, levels = LEV),
                      locus  = factor(locus,  levels = names(LOCI)))
pts <- pts %>% mutate(design = factor(design, levels = LEV),
                      locus  = factor(locus,  levels = names(LOCI)))

g <- ggplot(d, aes(n, power, colour = design)) +
  geom_hline(yintercept = 0.8, linewidth = 0.3, colour = "grey65", linetype = "22") +
  geom_line(aes(linetype = part, group = interaction(design, part)), linewidth = 0.8) +
  geom_point(data = pts, size = 2.2) +
  facet_wrap(~locus, ncol = 2) +
  scale_x_log10(breaks = c(100, 300, 1000, 3000, 10000, 30000, 100000),
                labels = c("100","300","1,000","3,000","10,000","30,000","100,000")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = expansion(0.015)) +
  scale_colour_manual(values = COL, name = NULL) +
  scale_linetype_manual(values = c(run = "solid", reach = "22"), guide = "none") +
  labs(x = "individuals/lines phenotyped", y = "power") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 8.5),
        legend.position = "bottom", legend.margin = margin(t = -4),
        legend.key.width = unit(20, "pt"))

ggsave(OUT, g, width = 7.2, height = 3.6, dpi = 300)
cat("\nWrote", OUT, "\n")
