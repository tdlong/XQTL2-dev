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
of phenotypic variance (Fig. 1b). Comparing non-overlapping 75 kb tiles against
the replicate error recovered from the odd/even split, the estimated contribution
exceeds that error at 83% of euchromatic tiles against 37% expected in its
absence, bounding the null fraction at 28% and implying that at least 72% of the
euchromatic genome carries non-zero heritability for lifespan (69–81% within a
single group). Individual contributions are nonetheless small. The chr3L locus
accounts for 4.6% of phenotypic variance in the group in which it is strongest
and 2.5% in the weakest — substantial for a longevity locus, but far from the
major QTL implied by a −log10 *P* of 207 — while the next largest regions reach
1.3–1.7% and the median euchromatic window 0.15% in males and 0.07% in females.
Resolving contributions of this magnitude required screening ~400,000 flies,
selecting the longest-lived 5.5% and sequencing 23,485 selected individuals;
effects in the second tier of Fig. 1a fall below the detection limit of designs
an order of magnitude smaller.

**Architecture differs between the sexes but is largely invariant to diet.** On
the X chromosome the median window contributes 0.43% of phenotypic variance in
males and 0.01% in females, males exceeding females at every window; on 2L males
exceed females at 79% of windows, whereas on 3R the direction reverses and
females exceed males at 71% (Fig. 1b). Hemizygosity accounts for part of the X
difference, exposing recessive variation that a second chromosome would mask, but
2L and 3R are ordinary autosomes. Diet behaves differently: the two sugar levels
alter lifespan substantially yet leave the distribution of contributing loci
nearly unchanged. Partitioning per-window heritability into a component shared by
all four groups and components differing by sex or by diet (Fig. 1c), 10.2% of
the genetic variance is sex-specific and 1.9% diet-specific, with no detectable
interaction; on autosomes alone, where hemizygosity cannot contribute, the
sex-specific fraction is 7.3%. Within the fifth of the genome carrying most of the
signal, where the ratio is best determined, sex and diet together account for
9.4% (95% CI 5.9–15.9%). This fraction is invisible to any analysis that pools
sexes and treats environment as a covariate — the standard construction wherever
environment cannot be manipulated — and it is not distributed uniformly: it is
concentrated at loci acting in one sex and not the other, which such an analysis
reports at a fraction of their true effect.

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

| quantity | value |
|---|---|
| chr3L 9.30 Mb, corrected h² (SY20 M / SY10 M / SY20 F / SY10 F) | 4.60 / 4.07 / 2.93 / 2.50 |
| next largest regions: 3R 8.68 / 3L 21.64 / 2L 14.85 / 2L 10.70 / X 10.10 | 1.67 / 1.34 / 1.46 / 1.16 / 1.54 |
| median window h², males / females | 0.15 / 0.07 |
| median window h² by group (SY10 M / SY20 M / SY20 F / SY10 F) | 0.18 / 0.14 / 0.09 / 0.06 |
| windows above 2× / 5× / 10× / 20× the median | 35 / 17 / 4.0 / 0.5 % |

**Fraction of the genome with non-zero h²** — on 1,537 non-overlapping euchromatic
tiles (every 15th window), comparing the shared component against the replicate
error from the odd/even split. Both are 1 df against 4 df, so under the null the
difference is positive only 37.4% of the time, not 50%; the excess over that
bounds the null fraction.

| quantity | value |
|---|---|
| tiles where the shared component exceeds replicate error | 82.7% euchromatin, 83.0% genome-wide |
| expected under no heritability (simulation, 2×10⁶ draws) | 37.4% |
| implied upper bound on the null fraction | 28% |
| implied lower bound on tiles with h² > 0 | 72% |
| same within one group (SY20 M / SY10 M / SY20 F / SY10 F) | 80.5 / 78.5 / 75.5 / 68.6 % |

The bound assumes a window with real signal is detected with probability 1. If
detection is imperfect the null fraction is smaller still, so 72% is conservative.

**Design scale**

| quantity | value |
|---|---|
| cages with recorded selection counts | 40 of 48 |
| selected flies sequenced in those cages | 23,485 |
| implied flies screened | ~408,000 |
| selected proportion: median (range) | 5.5% (3.1–10.8) |

**Sex and diet (Fig. 1c)**

| quantity | value |
|---|---|
| median window h², M vs F: chrX | 0.429 vs 0.006, M higher at 100% of windows |
| 2L | 0.128 vs −0.015, M higher at 79% |
| 2R | 0.044 vs −0.011, M higher at 67% |
| 3L | 0.165 vs 0.146, M higher at 69% |
| 3R | 0.113 vs 0.222, M higher at 29% |
| split, all chromosomes (shared / sex / diet / sex×diet) | 88.5 / 10.2 / 1.9 / −0.6 % |
| split, autosomes | 90.9 / 7.3 / 2.2 / −0.5 % |
| split, X only | 59.4 / 44.3 / −1.8 / −1.9 % |
| top 20% of genome by Wald: sex+diet share | 9.4% [5.9, 15.9] |

The shares (10.2%, 7.3%, 9.4%) are ratios of sums over the same windows, so the
fifteenfold window overlap cancels from numerator and denominator. Absolute
totals do not survive it, which is why none are given.

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
