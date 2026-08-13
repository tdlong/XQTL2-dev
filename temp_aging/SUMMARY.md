# AGE_SY — partitioning heritability into shared, sex and diet components

Working notes for the aging dataset: what we are asking, how the numbers in
Figure 1 were produced, and every assumption they rest on.

## The question

A polygenic risk score is fitted after sex and environment have been regressed
out, so it estimates only the component of genetic effect that is common to
everyone. That procedure cannot distinguish "there is no sex- or
environment-specific component" from "there is one and we removed it." How big
is the part that gets discarded?

AGE_SY can answer that because it is stratified by construction: lifespan
selection applied to a synthetic population in a 2x2 of sex (female, male) and
diet (SY10, SY20), 12 replicate cages per cell. Environment is assigned rather
than observed, so there is nothing to correct for and no confounding with
ancestry or with anything correlated with geography — which is the reason the
same quantity is hard to estimate in human cohorts.

## The data

- 8 founders, ~75 kb haplotype windows sliding in 5 kb steps, smoothed at 100 kb.
- Four scans: `AGE_SY10_F`, `AGE_SY20_F`, `AGE_SY10_M`, `AGE_SY20_M`.
- Each scan contrasts the selected (long-lived) pool against its unselected
  control, and reports a Wald test of haplotype frequency change plus Falconer
  and Cutler heritability.
- **Cutler h2 is a percentage** — `Heritability()` multiplies by 200 = 2 x 100.
  A reported 0.7 is 0.7% of the variance at that window. Genome median 0.68%,
  chr3L peak ~4.9%. That is the expected scale per 75 kb window for a polygenic
  trait.

### Getting an error term: the odd/even split

The pipeline reports one h2 per scan, averaged over replicates, so there is no
error term for a partition (XQTL2 #33; closed without change, per-replicate
values are not wanted in the output). We manufactured one by running each of the
four treatments twice — once on odd replicates, once on even — reusing the
existing haplotypes. That gives **8 measures per window**: 4 treatments x 2
halves, with the four within-cell differences supplying a 4 df pure error.

Odd/even rather than R1-6 vs R7-12 because the replicates arrived over several
months; a block split would put batch into the error term.

Design files keep their real `REP` values (odd = 1,3,5,7,9,11). This is only
safe because XQTL2 #32 was fixed first — before it, `Heritability()` matched
replicates by position rather than label and those files would have silently
produced h2 from 3 of the 6 replicates with mismatched selection proportions,
while the Wald score looked perfectly healthy. **The first run of these scans hit
exactly that**: the cluster's pipeline checkout was behind #32, and every number
from that pass was discarded.

## The partition

At each window the 8 measures are decomposed as a balanced 2x2 with 2 replicates
per cell:

```
sum(y^2) = n*ybar^2 + SS_sex + SS_diet + SS_sex:diet + SS_rep
```

with `n*ybar^2 = 8*ybar^2` (1 df), each treatment term 1 df, and `SS_rep` on
4 df. The replicate term is pure error; since it carries 4 df and the others 1,
it is subtracted on **mean** squares — each term becomes `MS - SS_rep/4`.
Verified against `aov()` at chr3L 9,400,000, and the parts sum to `sum(y^2)`
to 9e-14.

Terms are also reported in h2 units as `sign * sqrt(|MS - MS_rep| / 8)`, since
every term satisfies `SS = 8 * deviation^2`. Panel C is in those units. They are
**not floored at zero**: a term below the replicate error goes negative, and
that is the signal that there is nothing at that window.

### The h2 floor, and why it is calibrated

h2 squares an estimated quantity, so noise inflates it: `E[h2] = true + b`, and
`b` does not shrink with replicates because it is a within-replicate estimation
bias identical in every replicate. This is why h2 is positive everywhere
(minimum 0.099 over 200k measures). XQTL2 #34 now reports `b` per window
(`Cutl_H2_bias`), computed by Gauss-Hermite quadrature of `E[Affect^2]` through
the penetrance clamp.

The reported `b` is well calibrated where it is small but badly overstated where
it is large. Measured on windows the Wald score calls null — where the
frequencies did not move, so true h2 ~ 0 and the observed h2 **is** the floor:

| raw b | h2/b at null windows |
|---|---|
| < 0.5 | 1.37 |
| 0.5–0.75 | 0.90 |
| 0.75–1 | 0.54 |
| 1–1.5 | 0.39 |
| 1.5–2 | 0.35 |
| > 2 | 0.16 |

Left uncorrected this drove corrected h2 to −6 at a few Mb, which is impossible.
We therefore fit an isotonic (monotone) regression of h2 on b over the
Wald < 2 windows and use the fitted value as the floor. The fitted curve is
nearly flat — raw b from 0.13 to 8.29 maps to 0.47–0.64 — so in practice the
correction is subtracting a near-constant ~0.55, and **the spatial structure in
the reported b is not real**.

Shrinkage toward a prior was considered and does not apply: the two halves agree
on b to 2.8%, giving an empirical-Bayes weight of 0.001. b is not noisily large
at the bad windows, it is precisely large, so it had to be recalibrated rather
than regularised.

## Assumptions

Everything above rests on these. They are listed in rough order of how much the
headline number depends on them.

1. **Low Wald implies true h2 ~ 0.** This defines the null used to calibrate the
   floor. It is the Wald test's own claim, but the calibration is fitted at
   Wald < 2 and applied at Wald > 100, which is an extrapolation. Mitigated by
   the fitted curve being flat, so there is little to extrapolate.
2. **The floor is estimable at all.** Without any floor removal the autosomal
   non-main share is 2.9%; with the calibrated floor it is 9.1%. The correction
   roughly triples the headline number. Range under +-25% error in the floor:
   6.8% to 11.4%.
3. **Odd and even halves are independent.** Required for `SS_rep` to be pure
   error. True for sampling and for the flies, since the cages are separate.
4. **b is common across the four cells.** Only then does it cancel from the sex
   and diet contrasts. It nearly does — the 8 measures agree on b to 14% — but
   selected fly counts differ (418–1087 in SY20_M against 314–701 in SY10_F),
   so a small part of any sex or diet term could be differential inflation.
5. **h2 aggregates linearly within a window.** The window bundles all causal
   variants into 8 founder effects, so per-variant sex effects of opposite sign
   partly cancel before being squared. The window-level sex share is therefore a
   **lower bound** on the average per-variant share.
6. **4 df is enough.** Every error term rests on 4 df per window. XQTL2 #33
   would give 44; it was declined.
7. **The diet contrast is a nudge.** SY20 put 20–25% more flies in the selected
   tail than SY10 (18 of 24 replicate-by-sex pairs, paired p = 0.014). Real, but
   mild — the diet component should be read as a floor on what environment can
   do, not an estimate of it.
8. **One trait, one species, eight founders, truncation selection.** Nothing
   here establishes that the magnitude generalises.

## Results

The question is how well the non-main share can be estimated, and where. Four
things had to be right before the number meant anything, each of which was got
wrong first:

- **windows overlap 15-fold** (75 kb windows stepping 5 kb), so summing raw
  windows counts each piece of genome fifteen times. Every 15th window tiles it
  once.
- **blocks must be genetic, not physical.** h2 is smeared over ~9 Mb in the
  low-recombination regions flanking the centromeres and ~1 Mb in euchromatin,
  so a fixed Mb block counts one region as ten. A 2 cM block here runs 0.38 to
  6.38 Mb depending on where it sits.
- **heterochromatin is included but collapsed** to one telomeric and one
  centromeric block per arm. It is 7% of the sequence but carries 20% of the h2;
  split at 2 cM it would fragment into dozens of blocks and dominate the
  ranking. Excluding it altogether was the earlier error -- it is real data.
- **the X is included.** This table is about how well the ratio can be
  estimated, not about what transfers to another species, and the X contributes
  both heritability and signal.

h2 is rescaled so the genome totals 50%, the known broad-sense heritability of
longevity. That is a real number, not a convenience; the ratio does not depend
on it in any case.

`partition_by_wald_rank.R` produces this. Blocks are given as
euchromatic+heterochromatic.

| top % genome | blocks | Mb | h2 (of 50) | main | sex+diet | % sex+diet | CI |
|---|---|---|---|---|---|---|---|
| 5 | 8+0 | 13.0 | 17.5 | 16.1 | 1.33 | 7.6 | [3.6, 13.2] |
| 10 | 13+2 | 26.2 | 33.1 | 30.7 | 2.35 | 7.1 | [4.3, 14.8] |
| 20 | 26+3 | 38.7 | 38.6 | 35.0 | 3.61 | 9.4 | [5.9, 15.9] |
| 30 | 40+3 | 54.1 | 41.9 | 37.7 | 4.12 | 9.8 | [6.0, 16.3] |
| 50 | 68+4 | 74.9 | 46.3 | 41.6 | 4.64 | 10.0 | [4.4, 15.8] |
| 100 | 136+7 | 124.9 | 50.0 | 45.0 | 4.78 | 9.6 | [-25.9, 20.9] |

**Thirteen blocks -- 26 Mb, a fifth of the genome -- hold 33 of the 50 points.**
Across the top 10-50% the estimate is **7-10%, with intervals around [4, 16]**.

Going wider gains little and eventually breaks: the last half of the genome adds
under 4 points of heritability and takes the interval to [-25.9, 20.9], because
where main is near zero the ratio has a pole in it.

### The interval

It integrates two sources of uncertainty in one bootstrap, rather than reporting
the floor as a separate sensitivity table:

1. **which blocks you sampled** -- resample 8 cM groups
2. **the h2 floor** -- the floor is an isotonic fit of h2 on the reported bias
   over Wald-null windows, so each replicate resamples those windows in blocks
   and refits, carrying its own floor

The floor is the larger of the two and is why the interval blows up once
low-h2 blocks are included.

### What the blocks are made of

Only 5 of the 8 autosomal arm ends exist as blocks, and 7 of 10 with the X: the
scan does not reach chr2R's centromere-proximal end or chr3R's telomere-proximal
end, and chr2R's other end is 20 kb, smaller than one tile. So the
heterochromatin here is **almost entirely chromosome 3's centromere** --
chr3L_het_end (3.0 Mb, 4.14 points) and chr3R_het_start (3.1 Mb, 4.64 points),
both reading ~0% sex+diet. Chromosome 2's centromere contributes 0.17 points.
Whether chromosome 2's would behave the same way is untestable here.

The X's heterochromatic blocks go the other way: chrX_het_end reads 49.9% and
chrX_het_start 80.8% non-main, on small absolute contributions.

### Other properties of the partition

- The **interaction is negative under every treatment of the floor** -- below
  the replicate error everywhere. No sex-by-diet interaction across this diet
  range.
- Male Wald scores on the X are *higher* than female (median 1.03 vs 0.30, 4.2%
  vs 0.0% of windows clearing Wald 6) despite males being hemizygous and so
  sampled at half the chromosomes. Lower coverage would cut male power, so the X
  signal is not a precision artifact. Dosage compensation equalises expression,
  not allele count, so it does not remove the dominance exposure that makes the
  X's sex component large.

### What this experiment can and cannot resolve

It resolves roughly **ten to twenty independent regions**. That is the
denominator for any statement about the partition, and it is why the interval is
[4, 16] rather than anything sharper. More replicates or finer windows will not
tighten it; more genome or a trait with more loci would.

It cannot count loci or recover an effect-size distribution. Peak-calling is
censored both ways: the ten peaks above h2 = 0.75 blank out 19% of the genome
under any exclusion rule, and the noise floor (median MS_rep 0.020) sits exactly
where a power law says the counts should be growing. The chr3L QTL is 1.1 Mb
wide at half-max and spans 4.4 cM, so no fixed exclusion window separates loci
there and in the pericentromeric blocks at the same time.

## Are the terms real? (the sign test, and its non-obvious null)

A term below zero means it is smaller than the replicate error. It is tempting
to read "below zero half the time" as noise, but that is the wrong baseline.
MS_term carries 1 df and MS_rep carries 4, and a 1 df chi-square is strongly
right-skewed (median 0.45 sigma^2 against mean sigma^2), so a term that is PURE
NOISE sits below zero **62.6%** of the time, not 50%. Simulated over 2e6 draws.

Observed, with 95% intervals from resampling 67 independent 2 Mb blocks:

| term | % below zero | 95% CI | reading |
|---|---|---|---|
| main | 16.7% | [12.2, 21.1] | real |
| sex | 29.9% | [22.1, 36.8] | real |
| diet | 47.5% | [40.4, 55.4] | real, weakest |
| sex:diet | 67.4% | [61.1, 73.0] | nothing -- CI contains the null |

Autosomes only: sex 35.6% [28.1, 44.0], diet 42.1% [34.0, 50.4]; both still
exclude 62.6.

So the stretches where main or sex dip below zero in Figure 1C are expected even
for a strongly real term, and diet spending about half its time below zero is
evidence FOR diet, not against it. The interaction sitting on the null is the
one term the test calls empty, which agrees with every other view of it.

## Open

- Panels A and B of Figure 1 currently use the split-half scans averaged over
  halves (6 replicates each), not the 12-replicate scans. Point `LIVE_SCAN_DIR`
  in `make_figure1.R` at `process/AGE_SY/Scans` to swap them in.
- Worth reporting back to XQTL2: the #34 floor gets the level about right and
  the spatial variation wrong. That will affect every experiment using it.
- The diet term at chr3L 8.5 Mb sits in a region where the raw b was ~2, so it
  is the least trustworthy of the diet signals; chr2L 14–15 is cleaner.

## Files

| | |
|---|---|
| `make_figure1.R` | builds `Figure1_plot.png` |
| `Figure1_legend.md` | figure legend |
| `numbers_for_text.txt` | the partition numbers |
| `../scripts_oneoffs/AGE_SY/splithalf/` | design files, scan submission, gather, partition |
| `../process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz` | the 8 measures per window |
| `../process/AGE_SY_splithalf/H2_varcomp_by_window.txt.gz` | per-window partition |
