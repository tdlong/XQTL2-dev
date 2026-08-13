# Figure 1

**An XQTL of lifespan.** (A) Wald test of haplotype frequency change between the
selected long-lived pool and its unselected control. (B) Heritability, *h*²,
estimated from haplotype frequency shifts, expressed as a percentage of
phenotypic variance attributable to the window. (C) *h*² partitioned into the
component shared by all four treatments (main), the component that differs
between the sexes (sex), and the component that differs between diets (diet).

Four treatments throughout: SY10 and SY20 diets in females and males, 12
replicate cages each; 75 kb haplotype windows in 5 kb steps, smoothed at 100 kb.
Chromosome arms are concatenated left to right (X, 2L, 2R, 3L, 3R); dotted lines
separate arms and grey bands mark heterochromatin (dm6 euchromatin boundaries,
Huynh et al. 2023 *PLoS Genet* 19:e1010439, Table S2). All three panels share one
horizontal axis. The vertical axes are broken so the chr3L locus does not set the
scale for the rest of the genome, at 45 in (A), 2.5 in (B) and 0.75 in (C); the
chr3L peak reaches −log10(*P*) = 207.

## Numbers

The chr3L peak at 9.30 Mb reaches −log10(*P*) = 207.5 in SY20 males, 133.6 in
SY10 males, 73.1 in SY20 females and 49.8 in SY10 females (A), with peak *h*² of
4.85, 4.42, 3.09 and 2.51 respectively (B). At 9.34 Mb, with the floor removed,
*h*² is 3.11 and partitions as main 3.11, sex 0.90, diet 0.20 — sex and diet
together account for 7.8% of the heritability at that locus.

Genome-wide, main accounts for 88.5% of the genetic variance, sex 10.2%, diet
1.9% and sex-by-diet −0.6%. On autosomes alone these are 90.9%, 7.3%, 2.2% and
−0.5%; on the X, 59.4% and 44.3% for main and sex. The X's large sex component
reflects hemizygosity in males exposing recessive variation that a diploid
autosome masks, a sex-chromosome property rather than evidence of
sex-differential effects.

The signal is concentrated: 13 blocks of 2 cM — 26 Mb, a fifth of the genome —
carry 33 of the 50 percentage points of trait heritability, and across the top
10–50% of the genome ranked by Wald the sex-plus-diet share is 7–10% (95% CI
[4, 16]).

## The partition in (C)

At each window the four treatment *h*² values are rewritten as

    h²(treatment) = main ± sex ± diet ± interaction

with **main** the average of the four, **sex** half the difference between the
female and male averages, **diet** half the difference between SY20 and SY10,
and **interaction** the extent to which the diet effect differs between the
sexes. The decomposition is exact — adding the components back with the
appropriate signs returns the four treatment values.

Each treatment was measured twice, on odd and on even replicates, so every window
carries eight values rather than four. The disagreement of a treatment with
itself is pure replicate error, and it is subtracted from every component. A line
below zero therefore means that component is smaller than the replicate noise,
i.e. nothing is there. Components are plotted in *h*² units and smoothed with a
500 kb rolling mean. The sex-by-diet interaction is not drawn: it lies below the
replicate error at 67% of windows and its 95th percentile is 0.08.

*h*² is upwardly biased because it squares a quantity estimated with error, so it
comes out positive even where nothing is happening, and the bias does not shrink
with replication. That floor is estimated from windows where the Wald test
detects no frequency change — where true *h*² is ~0, so the observed *h*² is the
floor itself — and removed before the partition in (C). Panels (A) and (B) show
the scans as estimated, uncorrected.
