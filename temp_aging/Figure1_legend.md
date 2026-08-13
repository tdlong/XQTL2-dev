# Figure 1

**Lifespan selection in a 2x2 of sex and diet, and the partition of heritability
into shared, sex and diet components.**

Chromosome arms are concatenated left to right (X, 2L, 2R, 3L, 3R); dotted
vertical lines separate arms and grey bands mark pericentromeric and telomeric
heterochromatin (dm6 euchromatin boundaries from Huynh et al. 2023 *PLoS Genet*
19:e1010439, Table S2). All three panels share the same horizontal axis.

**(A)** Wald test of haplotype frequency change between the selected long-lived
pool and its unselected control, as -log10(*P*), for each of the four
treatments: SY10 and SY20 diets in females and males. One curve per treatment,
75 kb haplotype windows in 5 kb steps, smoothed at 100 kb.

**(B)** Cutler broad-sense heritability *h*² at the same windows, one curve per
treatment, same colours as (A). *h*² is expressed as a percentage of phenotypic
variance attributable to the window. Values are as estimated, with no correction
for the upward bias described below.

**(C)** *h*² partitioned into the component shared by all four treatments
(main), the component that differs between the sexes (sex), and the component
that differs between diets (diet). At each window the four treatments were
measured twice — on odd and on even replicates — giving 8 values, decomposed as
a balanced 2x2 with 2 replicates per cell. Each term is shown as
sign x sqrt((MS - MS_rep)/8), i.e. returned to *h*² units and directly
comparable with (B), after subtracting both the replicate mean square (pure
error, 4 df) and the estimation floor. Terms are not floored at zero: a value
below zero means the term is smaller than the replicate error, i.e. nothing is
there. Curves are smoothed with a 500 kb rolling mean. The vertical axis is
clipped at 1; main and sex exceed it at the chr3L peak and main exceeds it over
chr3L 20-26 Mb and chr3R 1-9 Mb. The sex-by-diet interaction is not drawn: it
lies below the replicate error at 67% of windows and its 95th percentile is
0.08.

*h*² is upwardly biased because it squares an estimated quantity, and the bias
does not shrink with replication. In (C) it is removed using a floor calibrated
against windows where the Wald test detects no frequency change, and where true
*h*² is therefore ~0 and the observed *h*² is the floor itself.

Autosomes only, the shared component accounts for 90.9% of the genetic variance,
sex for 7.3%, diet for 2.2%, and the sex-by-diet interaction for -0.5%. The X is
excluded from that summary: hemizygosity in males exposes recessive variation
that a diploid autosome masks, so its large apparent sex component (44.3%) is a
sex-chromosome property rather than evidence of sex-differential effects.
