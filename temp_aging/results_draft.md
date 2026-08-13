# Results (draft)

**Selection for long life moves haplotype frequencies over much of the genome.**
Each of the four sex-by-diet groups was scanned separately, contrasting the
long-lived tail of each of 12 cages against an unselected sample from the same
cage (Fig. 1a). A Bonferroni threshold over the ~1,670 independent 75 kb windows
that tile the genome falls at −log10 *P* = 4.5, and 45% of euchromatic windows
exceed −log10 *P* = 6 in at least one group; within a single group the figure is
24–32%. Significance therefore does not usefully divide this genome into hits and
non-hits, because too much of it qualifies. One region is far above the rest: a
locus at chr3L 9.3 Mb reaches −log10 *P* = 207.5 in high-sugar males and is the
strongest signal in all four groups. What the scan is not is a featureless
elevation. The 40 contiguous euchromatic regions above threshold have a median
width of 530 kb (1.1 cM), with a quartile range of 180 kb to 1.7 Mb, and the
narrowest of them are as narrow as the 75 kb window and 100 kb smoothing can
render — the resolution limit here is the analysis, not the extent of linkage
disequilibrium in a population that has been recombining freely for many
generations. Structure on that scale is easier to reconcile with many discrete
loci of small effect than with a strictly infinitesimal architecture, though the
scan alone cannot separate the two.

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
| windows / independent tiles / Bonferroni 0.05 threshold | 24,988 / 1,666 / −log10 *P* = 4.52 |
| euchromatic windows above 4.52 / 6 / 8 / 10 | 59.3 / 45.4 / 35.6 / 28.7 % |
| above 6 within one group (SY10 M / SY10 F / SY20 F / SY20 M) | 31.6 / 26.0 / 26.2 / 24.1 % |
| max −log10 *P* per group (SY20 M / SY10 M / SY20 F / SY10 F) | 207.5 / 133.6 / 73.1 / 49.8 |
| contiguous euchromatic regions above −log10 *P* = 6 | 40 |
| their width: median / quartiles / max | 530 kb (1.09 cM) / 178 kb–1.7 Mb / 10.4 Mb |
| their peak height: median / upper quartile | 9.1 / 14.4 |

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
- **No count of loci.** Clumping the scan top-down with a ±1 cM exclusion returns
  102 "peaks", but six of the top seven are shoulders of the chr3L signal, which
  is 3.5 Mb wide above threshold. No fixed exclusion width works for both that
  locus and the sub-Mb regions elsewhere, so the number is an artifact of the
  choice. The 40 contiguous runs quoted above are a description of the scan, not
  a count of causal loci.
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
