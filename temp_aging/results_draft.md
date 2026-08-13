# Results (draft)

**Heritability for lifespan is spread across most of the genome.** Selecting the
longest-lived flies moved haplotype frequencies at windows throughout the genome
rather than at a handful of loci (Fig. 1a,b). After removing the estimation
floor, 78% of euchromatic windows still carry positive heritability, and the
shared component exceeds the replicate error at 83% of them. Almost all of these
windows are individually small: the median window accounts for 0.15% of
phenotypic variance in males and 0.07% in females, and only 4% of windows reach
ten times the median. One region stands apart — a locus at chr3L 9.3 Mb reaching
−log10 *P* = 207 and 4.6% of phenotypic variance in high-sugar males, thirty
times a typical window — but it has few companions: 0.5% of windows reach twenty
times the median. We do not report a genome-wide total. The 75 kb windows step 5
kb and so overlap fifteenfold, and each value is a point estimate carrying an
upward bias that does not shrink with replication, so summing them measures
neither. Everything below concerns the distribution of heritability across
windows and the way it splits, not its integral.

**Which parts of the genome matter depends on sex, and hardly at all on diet.**
Lifespan does not have one genetic basis shared by all flies. On the X
chromosome the median window carries 0.43% of phenotypic variance in males and
0.01% in females, and males are higher at every single window; on 2L males are
higher at 79% of windows, while on 3R the direction reverses and females are
higher at 71% (Fig. 1b). Part of the X difference is expected, since males carry
one copy and recessive variation that a second chromosome would mask is exposed
— but 2L and 3R are ordinary autosomes and that explanation does not reach them.
Diet behaves quite differently. The two sugar levels change lifespan
substantially, yet they barely change which parts of the genome contribute to
it. Splitting the heritability at each window into the part shared by all four
groups and the parts that differ by sex or by diet (Fig. 1c), 10.2% of the
genetic variance is specific to one sex and 1.9% to one diet, with no detectable
interaction between them; on the autosomes alone, where hemizygosity cannot
contribute, the sex share is 7.3%. Within the fifth of the genome carrying most
of the signal, where the split is best determined, sex and diet together account
for 9.4% (95% CI 5.9–15.9%). That is the part an analysis which pools the sexes
and treats environment as a covariate — the usual construction when environment
cannot be manipulated, as in human studies — is unable to see. It is not spread
thinly across the genome either: it sits at loci that act in one sex and not the
other, which such an analysis would report at a fraction of their real effect.

---

## Numbers behind the above

Per-window values are corrected heritability, in percentage points of phenotypic
variance, over euchromatic windows only.

| quantity | value |
|---|---|
| euchromatic windows with h² > 0 after floor removal | 77.7% (M 79.0, F 65.8) |
| windows where the shared component exceeds replicate error | 83.3% |
| median window h², males / females | 0.15 / 0.07 |
| median window h² by condition (SY10 M / SY20 M / SY10 F / SY20 F) | 0.18 / 0.14 / 0.06 / 0.09 |
| windows above 2× / 5× / 10× / 20× the median | 35 / 17 / 4.0 / 0.5 % |
| peak ÷ median within a condition | 21–33× |
| median window h², M vs F: chrX | 0.429 vs 0.006, M higher at 100% |
| 2L | 0.128 vs −0.015, M higher at 79% |
| 2R | 0.044 vs −0.011, M higher at 67% |
| 3L | 0.165 vs 0.146, M higher at 69% |
| 3R | 0.113 vs 0.222, M higher at 29% |
| chr3L peak Wald −log10 *P* (SY20 M / SY10 M / SY20 F / SY10 F) | 207.5 / 133.6 / 73.1 / 49.8 |
| chr3L peak corrected h² | 4.60 / 4.07 / 2.93 / 2.50 |
| split, all chromosomes (shared / sex / diet / sex×diet) | 88.5 / 10.2 / 1.9 / −0.6 % |
| split, autosomes | 90.9 / 7.3 / 2.2 / −0.5 % |
| split, X only | 59.4 / 44.3 / −1.8 / −1.9 % |
| top 20% of genome by Wald: sex+diet share | 9.4% [5.9, 15.9] |

The shares (10.2%, 7.3%, 9.4%) are ratios of sums over the same windows, so the
fifteenfold window overlap cancels from numerator and denominator and they are
not affected by it. Absolute totals are, which is why none are given.

## Changed from the earlier draft, and why

- **Per-condition heritability totals removed** (they read 70 / 57 / 41 / 32%).
  They were circular: window h² was summed, which overcounts fifteenfold, and
  then rescaled so the genome came to 50% — the number that was supposed to be
  the result. Nothing in this experiment estimates broad-sense heritability for
  lifespan overall. Replaced with per-window medians and the male:female
  comparison at matched windows, which need no integration.
- **"93% of H² lies outside the top ten windows" removed.** Same problem: a
  statement about a sum of overlapping point estimates. The concentration claim
  is now made by counting windows above multiples of the median.
- **Male:female contrasts restated per window.** The earlier chromosome ratios
  (66.9× on X, 9.0× on 2L) were ratios of sums, so they inherit whatever the
  denominator does when it approaches zero — the female median on 2L is in fact
  slightly negative, which is why that ratio was large. Median-versus-median plus
  the fraction of windows where males exceed females says the same thing without
  a pole in it.
- **"essentially the entire euchromatic genome" → "most".** 22% of euchromatic
  windows sit at or below zero after floor removal.
- **2R dropped from the female-biased claim.** It is mildly male-biased. 3R is
  female-biased and that part stands.
- **chr3L peak values corrected.** Not "~5 in one condition and >3 in the
  others": 4.60 and 4.07 in males against 2.93 and 2.50 in females.
- **Second paragraph rewritten for a reader who has not done this analysis.**
  The previous version led with the machinery of the decomposition. The point is
  that lifespan genetics differs between the sexes and barely differs between
  diets, and that a pooled analysis cannot recover the difference.

## Not written, because the analysis does not exist

The third paragraph of the earlier draft requires three things that have not been
done here:

1. a null model for the local heritability expected from polygenic aggregation,
   against which regions could be called as exceeding it;
2. fine-mapping of the chr3L peak — it is 1.1 Mb wide at half-maximum and spans
   4.4 cM, and nothing in this analysis supports a kb-scale interval or names a
   gene;
3. founder-effect decomposition at the peak — per-founder effects are computed
   inside `Heritability()` and discarded, so this needs the smoothed frequency
   files rather than the scan output.
