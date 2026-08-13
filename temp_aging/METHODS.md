# Methods

## Design

Lifespan selection was applied to a synthetic outbred population derived from 8
inbred founders, in a 2 x 2 of sex (female, male) and diet (SY10, SY20), with 12
replicate cages per treatment (R1–R12). Within each replicate the long-lived
tail was collected as the selected pool and contrasted against an unselected
control pool from the same cage. The fraction of the population falling in the
selected tail averaged 0.047 (SY10 female), 0.056 (SY20 female), 0.055 (SY10
male) and 0.069 (SY20 male); SY20 exceeded SY10 in 18 of 24 replicate-by-sex
pairs (paired *t* on the log ratio, *P* = 0.014).

Control samples are shared between diets within a sex — the SY10 and SY20 design
files for a given sex reference the same 12 `Con` pools — so the diet and
sex-by-diet contrasts have control sampling error in common, while the sex
contrast does not.

## Sequencing, haplotype inference and scans

Pooled sequencing, reference/alternate counting, haplotype frequency estimation
and the genome scan were performed with the XQTL2 pipeline
(https://github.com/tdlong/XQTL2). Founder haplotype frequencies were estimated
in 75 kb windows sliding in 5 kb steps and smoothed over 100 kb. Each of the four
treatments was scanned separately (`run_scan.sh`: `smooth_haps` → `hap_scan` →
`concat_scans`), contrasting selected against control pools across the 12
replicates, and yielding per window a Wald test of haplotype frequency change
(−log10 *P*), Falconer and Cutler broad-sense heritability, and genetic position
in cM.

Cutler *h*² is computed per replicate as 200 × Σ_f C_f · Affect_f², where C_f is
the control frequency of founder *f*, Affect_f = Φ⁻¹(1−P) − Φ⁻¹(1−Pen_f) with
Pen_f = P·Z_f/C_f clamped to [P/2, 2P], P the selected proportion and Z_f the
selected-pool frequency; the per-replicate values are then averaged. The factor
of 200 makes *h*² a percentage, so a reported 0.7 is 0.7% of phenotypic variance
attributable to that window.

Two pipeline defects were identified and fixed before the analysis reported here
(XQTL2 issues #32 and #34); a third was reported and declined (#33). All scans
used a pipeline build containing the #32 fix. The split-half scans additionally
contain the #34 bias columns; the 12-replicate scans used for Figure 1A,B predate
#34 and carry none, which does not affect those panels.

## Replicate error: the odd/even split

The pipeline reports one *h*² per scan, averaged over replicates, so a single
scan provides no error term for a variance partition. Each of the four
treatments was therefore scanned twice more, once on odd replicates
(R1,3,5,7,9,11) and once on even (R2,4,6,8,10,12), reusing the existing
haplotype estimates and changing only the design file. Odd/even rather than a
block split because replicates were collected over several months; a block split
would confound batch with the error term.

Design files retain their true `REP` labels. This is safe only because XQTL2 #32
had been fixed: previously `Heritability()` matched replicates to selection
proportions by position rather than by label, so a subsetted design silently
produced *h*² from 3 of 6 replicates with mismatched proportions while the Wald
score appeared normal.

This yields **8 measures per window**: 4 treatments × 2 halves.

## Partition

At each window the 8 values are decomposed as a balanced 2 × 2 with n = 2
replicates per cell, uncorrected for the mean:

```
Σy² = n·ȳ² + SS_sex + SS_diet + SS_sex:diet + SS_rep
```

with n·ȳ² = 8ȳ² (1 df), each treatment term 1 df, and SS_rep — the sum over
cells of (odd − even)²/2 — on 4 df. Terms were verified against `aov()`; the
components sum to Σy² to 9 × 10⁻¹⁴.

Because SS_rep carries 4 df and the others 1, the error is subtracted on **mean**
squares: each term becomes MS − SS_rep/4. Terms are not truncated at zero; a
negative term indicates a component smaller than the replicate error. For display
and for reporting in *h*² units, each term satisfies SS = 8·deviation², so
components are given as sign × √(|MS − MS_rep|/8), with the main term signed by
ȳ rather than by the component (8ȳ² is a square, so an over-corrected window
would otherwise present as a large positive shared effect).

## The *h*² floor

Cutler *h*² squares a quantity estimated with error, so E[*h*²] = true + b with
b > 0, and b does not diminish with replication: it is a within-replicate bias
identical in every replicate, so averaging reduces the variance of the estimate
but not the bias. Empirically *h*² never falls below 0.114 across the 199,904 measures
(24,988 windows × 4 treatments × 2 halves).

XQTL2 #34 reports b per window (`Cutl_H2_bias`), obtained by 7-point
Gauss–Hermite integration of E[Affect²] over Pen ~ N(Pen_raw, Var(Pen)) with the
penetrance clamp inside the integrand. The reported b is well calibrated where it
is small but overstated where it is large: measured on windows the Wald test
calls null, the ratio *h*²/b runs 0.90 for b in 0.5–0.75 and 0.16 for b > 2.
Subtracting it unmodified drives corrected *h*² to −6 at some windows, which is
inadmissible.

The floor was therefore **recalibrated against an internal null**. Windows with
Wald −log10 *P* < 2 have no detectable frequency change, so true *h*² ≈ 0 there
and the observed *h*² is the floor itself. An isotonic (monotone non-decreasing)
regression of *h*² on b was fitted over those windows and applied genome-wide.
The fitted relationship is nearly flat — raw b spanning 0.13 to 8.29 maps to
0.47–0.64 — so in practice this subtracts a near-constant ≈0.55 and implies that
the spatial structure in the reported b is not real. This calibration is ours,
not the pipeline's, and is the largest single source of uncertainty in the
results; it is propagated into all reported intervals (below).

Empirical-Bayes shrinkage of b was considered and rejected: the two halves agree
on b to 2.8%, giving a shrinkage weight of 0.001. b is not noisily large at the
problem windows, it is precisely large, and so required recalibration rather
than regularisation.

The floor cancels exactly from the sex, diet and sex-by-diet contrasts, which are
differences between cells; it enters only the main term.

## Blocks

Windows overlap 15-fold (75 kb sliding 5 kb), so summing *h*² across raw windows
counts each interval fifteen times. Every 15th window was taken, tiling the
genome once in 75 kb non-overlapping tiles.

Tiles were aggregated into **2 cM blocks**. Genetic rather than physical distance
because *h*² is smeared over ~9 Mb in the low-recombination regions flanking the
centromeres and ~1 Mb in euchromatin, so a fixed physical block treats one region
as many: a 2 cM block here spans 0.38–6.38 Mb depending on location. cM was
obtained as in the pipeline's `add_genetic()`: the *D. melanogaster* r6 genetic
map (`helpfiles/flymap.r6.txt`) smoothed with a normal kernel of 3 Mb bandwidth
and interpolated by spline.

Heterochromatin was **included but collapsed** to one telomere-proximal and one
centromere-proximal block per arm, using dm6 euchromatin boundaries (Huynh et al.
2023 *PLoS Genet* 19:e1010439, Table S2). Heterochromatin is 7% of the scanned
sequence but carries 20% of the *h*²; split at 2 cM it would fragment into many
blocks and dominate any ranking. Seven of the ten possible arm-end blocks exist:
the scan does not extend into chr2R's centromere-proximal or chr3R's
telomere-proximal heterochromatin, and chr2R's remaining end (20 kb) is smaller
than one tile. The heterochromatin represented is therefore predominantly the
chromosome 3 centromere.

The X chromosome is included.

## Scaling

Summed window *h*² is not a heritability: neighbouring windows measure
overlapping components of the same signal and the genome sums to several hundred
percent. Absolute columns were therefore rescaled so the genome totals 50%, the
broad-sense heritability of *Drosophila* longevity. This assumes the overcounting
is uniform across the genome, which is unlikely to hold exactly where haplotypes
are long. All ratios are invariant to the scaling.

## Intervals

Reported 95% intervals are percentile bootstraps over 300 replicates,
integrating two sources of uncertainty simultaneously:

1. **block sampling** — 8 cM groups resampled with replacement;
2. **the floor** — the Wald-null windows are resampled in 2 Mb blocks and the
   isotonic floor refitted, so each replicate carries its own floor.

The floor dominates, and is why intervals widen sharply once blocks with *h*²
near zero are admitted: the ratio has a pole where the main term passes through
zero.

## What the design does not resolve

The experiment resolves roughly 10–20 independent regions, and that — not the
number of windows or replicates — sets the width of all reported intervals.
Locus counting and effect-size distributions are not recoverable: peak calling is
censored in both directions, since the ten peaks with *h*² > 0.75 exclude 19% of
the genome under any separation rule, while the replicate noise floor (median
MS_rep = 0.021) lies exactly where small-effect counts would be expected to grow.
The chr3L peak is 1.1 Mb wide at half-maximum and spans 4.4 cM, so no single
exclusion width separates loci there and in the pericentromeric blocks
simultaneously.

## Code and data

Analysis code is in `temp_aging/` and `scripts_oneoffs/AGE_SY/splithalf/` of the
XQTL2-dev repository:

| | |
|---|---|
| `make_splithalf_designs.R` | the 8 odd/even design files |
| `run_splithalf_scans.sh` | submits the 8 scans, reusing existing haplotypes |
| `gather_splithalf_H2.R` | assembles the 8 measures per window |
| `varcomp_H2.R` | the per-window partition |
| `temp_aging/make_4scan_df.R` | collapses the four 12-replicate scans (run on cluster) |
| `temp_aging/make_figure1.R` | Figure 1 |
| `temp_aging/partition_by_wald_rank.R` | the block table and its intervals |

Analyses were run in R 4.0.3 with tidyverse; the pipeline runs under R 4.2.2.
