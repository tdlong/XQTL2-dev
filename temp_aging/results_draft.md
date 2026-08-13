# Results (draft)

**Selection for long life moves haplotype frequencies over much of the genome.**
Each of the four sex-by-diet groups was scanned separately, contrasting the
long-lived tail of each of 12 cages against an unselected sample from the same
cage (Fig. 1a). Windows overlap fifteenfold and are in linkage disequilibrium
besides, so the tests are far from independent and no multiple-testing correction
is attempted; taking −log10 *P* = 5 as a working cutoff, 55% of the euchromatic
genome exceeds it in at least one group, 29% exceeds 10 and 20% exceeds 15, and
within a single group the corresponding figures are 28–38%, 17–19% and 12–15%. A
threshold does not usefully divide this genome into hits and non-hits, because
too much of it qualifies. One region stands far apart: a locus at chr3L 9.3 Mb
reaching −log10 *P* = 207.5 in high-sugar males, the strongest signal in all four
groups. The better peaks also give a sense of the resolution available. This is
not a single-generation cross but an advanced-generation one, in which repeated
rounds of recombination have broken the founder haplotypes into short blocks.
Across the twelve well-separated peaks above −log10 *P* = 15, the interval within
2 units of the peak has a median width of 0.5 cM, and 0.11–0.84 cM for all but
one of them — corresponding to anywhere between 50 kb and 1.4 Mb of sequence,
depending on the local recombination rate. The genetic width of these intervals
is consistent; the physical width is not, and tracks recombination. The narrowest
of them are at or below the 75 kb window and 100 kb smoothing, so what limits
resolution here is the analysis rather than the extent of linkage disequilibrium
in the population.

**On a scale of trait variance, none of it is large.** The Wald test asks whether
a haplotype frequency moved, not how much of the trait the region controls, and
at this sample size a very small shift is very significant. To put the scan on a
variance scale we use the founder frequencies themselves. The population descends
from eight inbred founders, so at every window we estimate the frequency of each
founder haplotype in the selected flies and in the control from the same cage. A
founder haplotype over-represented among the long-lived is one whose carriers
were more likely to reach the selected tail, and the ratio of its two frequencies
is the penetrance of that haplotype for the selection criterion. Treating lifespan
as an underlying liability with a threshold at the selected proportion, each
penetrance becomes a displacement in mean liability for carriers of that
haplotype, and the variance in displacement among the eight founders — weighted
by their frequencies — is the share of phenotypic variance the window accounts
for (Methods). Figure 1b is that quantity. On it, the chr3L locus accounts for
4.6% of phenotypic variance in the group where it is strongest and 2.5% where it
is weakest: large for a longevity locus, but not the classical major QTL that a
−log10 *P* of 207 might suggest. The next largest regions reach only 1.3–1.7%,
and the median euchromatic window 0.15% in males and 0.07% in females. Effects of
that size are why the experiment is the size it is — roughly 400,000 flies
screened across the 40 cages with recorded counts, the longest-lived 5.5% taken,
and 23,500 selected flies sequenced. The second tier of Fig. 1a is visible
because of that scale rather than because the effects are big; in a conventional
experiment an order of magnitude smaller it would sit under the noise.

**Which parts of the genome matter depends on sex, and hardly at all on diet.**
Lifespan does not have one genetic basis shared by all flies. On the X chromosome
the median window carries 0.43% of phenotypic variance in males and 0.01% in
females, and males are higher at every single window; on 2L males are higher at
79% of windows, while on 3R the direction reverses and females are higher at 71%
(Fig. 1b). Part of the X difference is expected, since males carry one copy and
recessive variation that a second chromosome would mask is exposed — but 2L and
3R are ordinary autosomes and that explanation does not reach them. Diet behaves
quite differently. The two sugar levels change lifespan substantially, yet they
barely change which parts of the genome contribute to it. Splitting the
heritability at each window into the part shared by all four groups and the parts
that differ by sex or by diet (Fig. 1c), 10.2% of the genetic variance is
specific to one sex and 1.9% to one diet, with no detectable interaction between
them; on the autosomes alone, where hemizygosity cannot contribute, the sex share
is 7.3%. Within the fifth of the genome carrying most of the signal, where the
split is best determined, sex and diet together account for 9.4% (95% CI
5.9–15.9%). That is the part an analysis which pools the sexes and treats
environment as a covariate — the usual construction when environment cannot be
manipulated, as in human studies — is unable to see. It is not spread thinly
across the genome either: it sits at loci that act in one sex and not the other,
which such an analysis would report at a fraction of their real effect.

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
| euchromatic windows with h² > 0 after floor removal | 77.7% (M 79.0, F 65.8) |
| windows where the shared component exceeds replicate error | 83.3% |

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
