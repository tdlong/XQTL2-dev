# Methods

**Design.** Lifespan selection on a synthetic outbred population from 8 inbred
founders, in a 2 × 2 of sex × diet (SY10, SY20), 12 replicate cages per
treatment. In each cage the long-lived tail (mean 4.7–6.9% of the population)
was contrasted against an unselected control from the same cage. Controls are
shared between diets within a sex.

**Scans.** Haplotype frequencies were estimated with the XQTL2 pipeline
(github.com/tdlong/XQTL2) in 75 kb windows sliding 5 kb, smoothed over 100 kb.
Each treatment was scanned separately, giving per window a Wald test of
frequency change and Cutler broad-sense heritability *h*² = 200 Σ_f C_f·Affect_f²,
where Affect_f = Φ⁻¹(1−P) − Φ⁻¹(1−P·Z_f/C_f) with the penetrance clamped to
[P/2, 2P]; the factor 200 expresses *h*² as a percentage of phenotypic variance.

**Replicate error.** The pipeline averages *h*² over replicates, leaving no error
term, so each treatment was scanned twice more — odd and even replicates
separately, reusing the same haplotype estimates. This gives 8 values per window
(4 treatments × 2 halves). Odd/even rather than a block split because replicates
accrued over months.

**Partition.** At each window the 8 values were decomposed as a balanced 2 × 2
with 2 replicates per cell, uncorrected for the mean:
Σy² = 8ȳ² + SS_sex + SS_diet + SS_sex:diet + SS_rep, the first four on 1 df and
SS_rep (within-cell odd−even) on 4. Each term is reported as MS − SS_rep/4, and
in *h*² units as sign·√(|MS − MS_rep|/8). Terms are not truncated at zero.
Verified against `aov()`.

**Heritability floor.** *h*² squares an estimate, so E[*h*²] = true + b with b > 0
and independent of replication (minimum observed *h*² 0.114). XQTL2 reports b,
but it is overstated where large and subtracting it unmodified yields *h*² of −6.
It was recalibrated by isotonic regression of *h*² on b over windows with Wald
−log10 *P* < 2, where true *h*² ≈ 0 so the observed value is the floor; the fit
is near-flat (b of 0.13–8.29 maps to 0.47–0.64). b cancels from the sex, diet and
interaction contrasts and enters only the main term.

**Blocks.** Windows overlap 15-fold, so every 15th was taken to tile the genome
once, then aggregated into 2 cM blocks — genetic distance because *h*² spreads
over ~9 Mb pericentromerically and ~1 Mb in euchromatin. Heterochromatin (dm6
boundaries, Huynh et al. 2023) was collapsed to one telomeric and one centromeric
block per arm. Absolute values were rescaled so the genome totals 50%, the
broad-sense heritability of *Drosophila* longevity; ratios are unaffected.

**Intervals.** Percentile bootstrap, 300 replicates, resampling 8 cM block groups
and refitting the floor on resampled null windows, so both sources enter together.

**Code.** `temp_aging/` and `scripts_oneoffs/AGE_SY/splithalf/` of XQTL2-dev.
