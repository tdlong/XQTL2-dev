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

**Heritability floor.** *h*² sums squared founder effects estimated with error,
so squaring inflates it and a window at which no frequency changed still returns
a positive value. The pipeline reports this inflation per measure as *b*,
obtained by integrating E[Affect²] over the sampling distribution of the
frequency estimates. *b* is not used directly: it is well behaved where small but
overstated where large, and subtracting it as reported drives corrected *h*² to
−6 at some windows.

The floor actually subtracted was obtained as follows. Windows were classified as
null when the largest Wald −log10 *P* across the four treatments was < 2; at such
windows the frequencies did not detectably move, so true *h*² ≈ 0 and the
observed *h*² estimates the floor. Over the null windows, *h*² and *b* were each
averaged across the 8 measures, and an isotonic (monotone non-decreasing)
regression of mean *h*² on mean *b* was fitted, giving F. Every one of the 8
measures at every window was then corrected as *h*² − F(*b*), using that measure's
own *b*, and the corrected values are what enter the partition. F is near-flat:
reported *b* spanning 0.13–8.29 maps to 0.47–0.64.

F enters only the main term; it cancels from the sex, diet and interaction
contrasts, which are differences between cells.

**Blocks.** Windows overlap 15-fold, so every 15th was taken to tile the genome
once, then aggregated into 2 cM blocks — genetic distance because *h*² spreads
over ~9 Mb pericentromerically and ~1 Mb in euchromatin. Heterochromatin (dm6
boundaries, Huynh et al. 2023) was collapsed to one telomeric and one centromeric
block per arm. Absolute values were rescaled so the genome totals 50%, the
broad-sense heritability of *Drosophila* longevity; ratios are unaffected.

**Intervals.** Percentile bootstrap, 300 replicates, resampling 8 cM block groups
and refitting the floor on resampled null windows, so both sources enter together.

**Code.** `temp_aging/` and `scripts_oneoffs/AGE_SY/splithalf/` of XQTL2-dev.
