# Results (draft)

**Selection for longevity shifts haplotype frequencies across much of the
genome.** Each sex × diet group was scanned independently, contrasting the
long-lived tail of each of 12 cages with an unselected sample from the same cage
(Fig. 1a). Windows overlap fifteenfold and are further correlated through linkage
disequilibrium, so no multiple-testing correction is meaningful; at a working
threshold of −log10 *P* = 5, 55% of the euchromatic genome responds in at least
one group, 29% exceeds 10 and 20% exceeds 15 (28–38%, 17–19% and 12–15% within a
single group). Thresholding does not partition this genome into hits and
non-hits. One region is exceptional: a locus at chr3L 9.3 Mb reaching −log10 *P*
= 207.5 in high-sugar males, the strongest signal in all four groups. Because the
mapping population is an advanced-generation rather than a single-generation
cross, founder haplotypes are broken into short blocks and the stronger peaks are
correspondingly narrow. Across the twelve well-separated peaks above −log10 *P* =
15, the interval within 2 units of the peak spans a median of 0.48 cM (0.11–0.84
cM for eleven of the twelve) but 50 kb to 1.4 Mb of sequence, physical width
tracking local recombination rate. The narrowest intervals approach the 75 kb
window and 100 kb smoothing, indicating that resolution is limited by the
analysis rather than by linkage disequilibrium in the population.

**Longevity is highly polygenic.** Frequency change establishes that a region
responded to selection but not the magnitude of its contribution, and at this
sample size very small shifts are highly significant. Because the population
derives from eight inbred founders and the selected proportion of each cage is
known, the squared changes in founder frequency at a window can be converted to
the heritability that window accounts for (Methods), placing the scan on a scale
of phenotypic variance (Fig. 1b). Within a character the two are driven by the
same frequency shifts and trace a single curve, so the significance threshold
transfers directly onto the heritability scale: fitting *h*² on Wald per
character and evaluating at −log10 *P* = 5 returns 0.72, 0.75, 0.77 and 0.83% of
phenotypic variance across the four groups, against ~0.5% where the Wald test
detects nothing. Windows contributing more than ~0.8% of phenotypic variance are
therefore significant, and 32–45% of the euchromatic genome exceeds that value in
each group; because −log10 *P* = 5 is a conservative cutoff rather than a
detection limit, a considerably larger fraction is suggestive. Individual
contributions are nonetheless small. The chr3L locus reaches 4.85% of phenotypic
variance in high-sugar males and 2.51% in low-sugar females — substantial for a
longevity locus, but far from the major QTL implied by a −log10 *P* of 207 — the
next largest regions reach 1.5–2.2%, and the great majority of the significant
genome lies just above 0.8%. Resolving contributions of this magnitude required
screening ~400,000 flies, selecting the longest-lived 5.5% and sequencing 23,485
selected individuals; effects in this range fall below the detection limit of
designs an order of magnitude smaller. This third of the genome is also where the
partition in Fig. 1c can be read, a ratio of components being meaningful only
where there is heritability to divide.

**Architecture differs between the sexes but is largely invariant to diet.**
Restricting to the significant genome — the 843 euchromatic tiles, 55% of the
total, exceeding −log10 *P* = 5 in at least one group — heritability is sharply
sex-dependent and its direction varies by chromosome arm. Males exceed females at
100% of significant tiles on the X, 82% on 2L and 81% on 3L, with median *h*²
ratios of 31, 4.4 and 1.7; on 3R the direction reverses and females exceed males
at 81% of tiles (ratio 0.54), while 2R shows no consistent difference. Hemizygosity
accounts for the X, exposing recessive variation that a second chromosome would
mask, but 2L, 3L and 3R are ordinary autosomes and their sex dependence, in both
directions, is not attributable to it. Diet behaves differently: the two sugar
levels alter lifespan substantially yet leave the distribution of contributing
loci nearly unchanged. Partitioning heritability across these tiles into a
component shared by all four groups and components differing by sex or by diet
(Fig. 1c), 10.9% of the genetic variance is sex-specific and 2.0% diet-specific,
with no detectable interaction; on autosomes alone the sex-specific fraction is
9.2%, and over the fifth of the genome carrying most of the signal sex and diet
together account for 9.4% (95% CI 5.9–15.9%). This fraction is invisible to any
analysis that pools sexes and treats environment as a covariate — the standard
construction wherever environment cannot be manipulated — and it is not
distributed uniformly: it is concentrated at loci acting in one sex and not the
other, which such an analysis reports at a fraction of their true effect.

---

## Numbers behind the above

Heritabilities are percentage points of phenotypic variance, floor-corrected,
over euchromatic windows unless noted.

**Scan and resolution (Fig. 1a, full 12-replicate scans)**

| quantity | value |
|---|---|
| euchromatic genome above −log10 *P* 5 / 10 / 15 / 20, any group | 54.7 / 28.7 / 20.3 / 16.7 % |
| same within one group, above 5 | 28.1–38.4% (SY20 M lowest, SY10 M highest) |
| same within one group, above 10 / 15 | 16.6–19.4 / 12.4–15.0 % |
| max −log10 *P* per group (SY20 M / SY10 M / SY20 F / SY10 F) | 207.5 / 133.6 / 73.1 / 49.8 |
| well-separated euchromatic peaks above 15 (±5 cM exclusion) | 12 |

2-unit support intervals for those twelve, sorted by peak height:

| peak | −log10 *P* | 2-unit drop, kb | 2-unit drop, cM |
|---|---|---|---|
| chr3L 9.30 | 207.5 | 50 | 0.20 |
| chr3R 8.66 | 42.1 | 855 | 0.31 |
| chr3L 21.64 | 38.2 | 1370 | 0.14 |
| chr2L 10.70 | 30.5 | 120 | 0.37 |
| chr2R 25.10 | 28.7 | 90 | 0.11 |
| chr3L 7.50 | 27.0 | 555 | 2.23 |
| chr2L 14.85 | 23.9 | 225 | 0.31 |
| chr2R 13.88 | 21.0 | 270 | 0.80 |
| chr2R 16.52 | 17.3 | 165 | 0.58 |
| chr3R 22.41 | 16.8 | 195 | 0.84 |
| chr3R 21.09 | 16.7 | 240 | 0.77 |
| chrX 10.09 | 15.9 | 175 | 0.65 |

Median 210 kb / 0.48 cM. The chr3L 7.50 entry sits exactly 5.0 cM from the 9.30
peak, at the edge of the exclusion, and may be part of the same signal — it is
also the one interval wider than 0.84 cM. The 50 kb interval at chr3L 9.30 is
narrower than the 75 kb window, so it should be read as "at the resolution
limit", not as a 50 kb localisation.

**Magnitude (Fig. 1b)**

**What value of h² is significant** (`h2_threshold.R`). Within a character, Wald
and *h*² trace a single curve, so the −log10 *P* = 5 cutoff transfers onto the
*h*² scale. Fit is a spline of *h*² on log(Wald) per character over 1,537
non-overlapping euchromatic windows; values are raw, as plotted in Fig. 1b.

| character | h² at Wald 2 | **h² at Wald 5** | se | h² at Wald 10 | % euchr. above h²(Wald 5) | % euchr. Wald > 5 |
|---|---|---|---|---|---|---|
| SY10 female | 0.58 | **0.72** | 0.007 | 0.89 | 32.8 | 28.3 |
| SY20 female | 0.59 | **0.75** | 0.007 | 0.89 | 31.9 | 31.3 |
| SY10 male | 0.68 | **0.77** | 0.010 | 1.09 | 45.0 | 38.8 |
| SY20 male | 0.63 | **0.83** | 0.012 | 1.03 | 34.8 | 27.9 |

The last two columns are two routes to the same cut and agree to within a few
points, as they must if the curve is monotone. Where they differ it is scatter
about the spline: the male *h*² distributions are wider, so more windows clear the
fitted *h*² than clear the Wald cut itself.

At windows where the character's own Wald is below 2 — the floor — raw *h*² is
0.52 / 0.50 / 0.56 / 0.55 (median) in the four groups.

**Effect sizes on the same raw scale**

| quantity | value |
|---|---|
| chr3L 9.30 Mb (SY20 M / SY10 M / SY20 F / SY10 F) | 4.85 / 4.42 / 3.09 / 2.51 |
| next largest regions: 3R 8.66 / X 10.09 / 2L 14.85 / 3L 21.64 / 2L 10.70 | 2.23 / 2.08 / 1.99 / 1.89 / 1.69 |
| floor-corrected equivalents at chr3L 9.30 | 4.60 / 4.07 / 2.93 / 2.50 |
| median floor-corrected window h², males / females | 0.15 / 0.07 |

**Design scale**

| quantity | value |
|---|---|
| cages with recorded selection counts | 40 of 48 |
| selected flies sequenced in those cages | 23,485 |
| implied flies screened | ~408,000 |
| selected proportion: median (range) | 5.5% (3.1–10.8) |

**Sex and diet (Fig. 1c)**

Restricted to the 843 significant euchromatic tiles (Wald > 5 in at least one
group), from `significant_regions.R`.

| arm | tiles | median h² M | median h² F | M/F | M > F at |
|---|---|---|---|---|---|
| chrX | 62 | 0.67 | 0.02 | 30.8 | 100% |
| chr2L | 169 | 0.22 | 0.05 | 4.44 | 82% |
| chr3L | 186 | 0.87 | 0.51 | 1.72 | 81% |
| chr2R | 172 | 0.16 | 0.18 | 0.87 | 59% |
| chr3R | 254 | 0.19 | 0.34 | 0.54 | 19% |

The "M > F at" column is **identical** with the floor removed and with it left in
(82 / 59 / 81 / 19 / 100 either way). The floor is a near-constant added to every
window; it cancels from a within-window comparison between groups, and where h²
is well clear of it the comparison does not depend on how it was estimated.

| split over significant tiles | main | sex | diet | sex×diet |
|---|---|---|---|---|
| all chromosomes | 87.5 | 10.9 | 2.0 | −0.4 |
| autosomes | 89.0 | 9.2 | 2.1 | −0.4 |
| top 20% of genome by Wald: sex+diet | | 9.4% [5.9, 15.9] | | |

The split, unlike the male:female comparison, does move if the floor is left in
(sex 3.8% rather than 10.9%), because the floor adds to every group equally and
so inflates the shared component. That is why Fig. 1c is floor-corrected and
Fig. 1b is not.

The shares are ratios of sums over the same windows, so the fifteenfold window
overlap cancels from numerator and denominator. Absolute totals do not survive
it, which is why none are given.

## What is deliberately not claimed

- **No genome-wide total heritability.** Windows are 75 kb stepping 5 kb, so they
  overlap fifteenfold, and each value is a point estimate with an upward bias
  that does not shrink with replication. Summing them measures neither. An
  earlier draft reported 70 / 57 / 41 / 32% for the four groups; those came from
  summing the overlap and then rescaling the genome to 50%, which is the number
  the exercise was supposed to produce.
- **No count of loci, and no multiple-testing correction.** Clumping the scan
  top-down with a ±1 cM exclusion returns 102 "peaks", but six of the top seven
  are shoulders of the chr3L signal. No fixed exclusion width works for both that
  locus and the sub-Mb regions elsewhere, so any such number is an artifact of the
  choice. The twelve peaks quoted above are chosen with a deliberately generous
  ±5 cM exclusion in order to have well-separated examples to measure widths on;
  they are not a claim that there are twelve loci. Likewise, an earlier draft
  quoted a Bonferroni threshold over ~1,670 tiles, which is meaningless when the
  windows overlap fifteenfold and are in LD on top of that. The cutoff is stated
  as a choice and a range of thresholds is given instead.
- **Male:female stated per window, not as a ratio of sums.** An earlier draft gave
  67× on X and 9× on 2L; those divide by a female sum near zero — 2L's female
  median is in fact slightly negative — so they are unstable. Median-versus-median
  plus the fraction of windows where males exceed females says the same thing
  without a pole in it.

## Not written, because the analysis does not exist

The third paragraph of the earlier draft requires three things that have not been
done here:

1. a null model for the local heritability expected from polygenic aggregation,
   against which regions could be called as exceeding it;
2. fine-mapping of the chr3L peak — it is 750 kb wide at half-maximum, 3.5 Mb
   above threshold, and nothing here supports a kb-scale interval or names a gene;
3. founder-effect decomposition at the peak — per-founder effects are computed
   inside `Heritability()` and discarded, so this needs the smoothed frequency
   files rather than the scan output.
