# Figure 2

**Seven peaks, and what the founder haplotypes do under them.** Each of the seven
labelled cells covers 1 Mb centred on a peak. The upper third is the Wald test in
all four treatments; the lower two thirds is the change in founder haplotype
frequency, selected minus control, for the treatment with the strongest signal at
that peak. The eighth cell asks whether these are longevity effects or founders
on their way out of the cages: at each peak we take the haplotype most enriched
in the long-lived pool (**most protective**) and the one most depleted (**most
susceptible**), and plot its frequency in the *control* pools against replicate.
Replicates accrued over months, so replicate order is time order.

Founder colours are fixed across the figure (XQTL2.Xplore `get_palette`);
founders averaging below 2.5% in the controls are drawn faded, since their
changes are noise. Each cell's locus label is drawn in that peak's colour, which
is also its colour in the eighth cell — no separate key. Axes are shared: Wald
0–30 and Δ frequency −0.06 to 0.06 everywhere except chr3L, which is on 0–200 and
−0.10 to 0.05 in both of its sub-panels. Tick labels appear only where a panel
introduces a scale.

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

Frequency changes are at the peak window, averaged over the 12 replicates. B6 is
the protective haplotype at four of the seven peaks; B1 is the susceptible one at
two.

## The control panel

The concern it addresses is that a founder haplotype could be leaving the cages
for reasons unrelated to lifespan, and would then appear depleted in any selected
pool. If so, the susceptible traces would slope down over replicates.

They do not. Regressing control frequency on replicate gives slopes between
−0.008 and +0.011 per replicate, and six of the fourteen reach *P* < 0.05 — but
that test treats a frequency trajectory as twelve independent points. Against a
random walk with the same per-step variance, **none of the fourteen reaches
*P* < 0.05**; chr2L:14.85 goes from *P* = 2 × 10⁻⁶ to *P* = 0.67. The protective
series average +0.0044 per replicate and the susceptible −0.0019, and paired
across the seven loci that difference gives *P* = 0.16.

There is no evidence of anything happening across replicates beyond drift, which
is what the panel is there to establish. The excursion at replicates 8–9, visible
in both halves at the 2L loci, moves several series together and looks like a
batch rather than a trend.

*Code: `temp_aging/make_figure2.R`, `make_zoom_means.R`, `control_drift.R`.*
