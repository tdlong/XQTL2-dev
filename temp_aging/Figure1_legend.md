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
replicate error at 67% of windows.

*h*² is upwardly biased because it squares a quantity estimated with error, so it
comes out positive even where nothing is happening, and the bias does not shrink
with replication. That floor is estimated from windows where the Wald test
detects no frequency change — where true *h*² is ~0, so the observed *h*² is the
floor itself — and removed before the partition in (C). Panels (A) and (B) show
the scans as estimated, uncorrected.
