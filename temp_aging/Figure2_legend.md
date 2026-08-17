# Figure 2

**Seven loci in detail.** Each labelled cell spans 1 Mb centred on a peak. Upper
subpanel: Wald test in all four treatments. Lower subpanel: change in founder
haplotype frequency, selected minus control, for the treatment with the strongest
signal at that peak. Eighth cell: frequency in the unselected controls across
replicates of the haplotype most enriched in the long-lived pool (most
protective) and the one most depleted (most susceptible), one line per locus.

Founder colours are fixed throughout; founders below 2.5% frequency in the
controls are faded. Each locus label is drawn in that peak's colour, which is its
colour in the eighth cell. Replicate order is time order. Axes are shared — Wald
0–30, Δ frequency −0.06 to 0.06 — except chr3L, at 0–200 and −0.10 to 0.05.

## Numbers

| peak | max −log10 *P* | treatment drawn | most protective | most susceptible |
|---|---|---|---|---|
| chrX:10.09 | 15.9 | SY10 male | B7 (+0.046) | B6 (−0.034) |
| chr2L:10.71 | 30.5 | SY20 male | B6 (+0.054) | B1 (−0.031) |
| chr2L:11.95 | 28.9 | SY10 male | B6 (+0.061) | B2 (−0.022) |
| chr2L:14.85 | 23.9 | SY10 male | B6 (+0.027) | AB8 (−0.048) |
| chr2R:13.88 | 21.0 | SY20 female | B3 (+0.043) | B5 (−0.035) |
| chr2R:16.52 | 17.3 | SY20 female | B2 (+0.023) | B1 (−0.036) |
| chr3L:9.31 | 207.5 | SY20 male | B6 (+0.038) | B3 (−0.109) |

Changes are at the peak window, averaged over 12 replicates. B6 is protective at
four of seven peaks and susceptible at chrX:10.09.

## The eighth cell

Tests whether a haplotype is simply leaving the cages, which would make it look
depleted in any selected pool. If so, the susceptible traces would slope down.

They do not. Six of fourteen slopes reach *P* < 0.05, but that treats a frequency
trajectory as twelve independent points. Against a random walk of the same
per-step variance, none does — chr2L:14.85 goes from *P* = 2 × 10⁻⁶ to 0.67.
Protective series average +0.0044 per replicate, susceptible −0.0019; paired
across loci, *P* = 0.16. Nothing departs from drift. The excursion at replicates
8–9 moves several series together and looks like a batch.

*Code: `temp_aging/make_figure2.R`, `make_zoom_means.R`, `control_drift.R`.*
