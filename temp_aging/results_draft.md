# Results (draft)

**Longevity architecture is continuously polygenic.** Broad-sense heritability for longevity is elevated above the noise floor across much of the euchromatic genome in every sex–diet combination (Fig. 1b). After removing the estimation floor, 78% of euchromatic windows retain positive H², and the main effect exceeds the replicate error at 83% of windows. Integrated over the mapped genome, total broad-sense H² is 70% in low-sugar males, 57% in high-sugar males, 41% in low-sugar females and 32% in high-sugar females — approximately twofold higher in males than females in both diets. In every condition the majority of this variance comes from a continuously elevated polygenic background rather than a few large-effect loci: 93% of total H² lies outside the ten highest-heritability windows and 80% outside the top fifty. The polygenic background is therefore a spatially resolved object rather than a threshold artifact, and it accounts for the bulk of the genetic variance contributing to longevity in this population.

**The polygenic background reorganizes across sex and diet.** The background is not a single shared architecture uniformly rescaled across contexts; distinct chromosomal regions contribute preferentially to male or female heritability (Fig. 1b). The X chromosome and 2L are strongly male-biased along their lengths, with integrated H² 67-fold and 9-fold higher in males than females respectively, while 3R is female-biased (male:female 0.67). Male bias on the X is expected in part from hemizygous exposure of partially recessive variation, but the parallel bias on 2L cannot be attributed to that mechanism and points to sex-differential contribution to longevity variance from an autosomal arm. Partitioning H² at each window into a component shared by all four treatments and components differing by sex or by diet (Fig. 1c) shows that 10.2% of the genetic variance is sex-specific and 1.9% diet-specific, with no detectable sex-by-diet interaction. Restricted to the fifth of the genome carrying two-thirds of the heritability, where the ratio is best estimated, sex and diet together account for 7–10% (95% CI 4–16%). A marginal analysis that regresses out sex and environment — the standard construction in human GWAS, where factorial manipulation of environment is not possible — estimates only the shared component and is structurally blind to this fraction.

---

## Numbers behind the above

| quantity | value |
|---|---|
| euchromatic windows with H² > 0 after floor removal | 77.7% |
| windows where main exceeds the replicate error | 83.3% |
| median corrected H² per euchromatic window | 0.136 |
| integrated H², SY10 male / SY20 male / SY10 female / SY20 female | 70 / 57 / 41 / 32 % |
| H² outside the top 10 / 50 / 100 windows | 93 / 80 / 66 % |
| integrated H² male:female — chrX / 2L / 2R / 3L / 3R | 66.9 / 9.0 / 1.38 / 1.68 / 0.67 |
| partition, all chromosomes (main / sex / diet / sex×diet) | 88.5 / 10.2 / 1.9 / −0.6 % |
| partition, autosomes | 90.9 / 7.3 / 2.2 / −0.5 % |
| partition, X only | 59.4 / 44.3 / −1.8 / −1.9 % |
| top 20% of genome by Wald: share of h², sex+diet | 66% of h², 9.4% [5.9, 15.9] |
| chr3L peak H² (SY20 M / SY10 M / SY20 F / SY10 F) | 4.60 / 4.07 / 2.93 / 2.50 |

Integrated H² percentages carry the scaling described in Methods — summed window
H² overcounts, so absolute values are rescaled to a genome total of 50%. The
male:female ratios are unaffected by that scaling; the absolute percentages
inherit it.

## Changed from the earlier draft, and why

- **"essentially the entire euchromatic genome" → "much of".** After floor
  removal 22% of euchromatic windows sit at or below zero. 78% is a substantial
  polygenic background but not the whole genome, and "significantly elevated"
  would require a per-window test that has not been run.
- **2R dropped from the female-biased claim.** It is mildly male-biased
  (male:female 1.38). 3R is female-biased (0.67) and that part stands. ("The left
  arm of 3R" also does not parse — 3R is an arm.)
- **chr3L peak values corrected.** Not "~5 in one condition and >3 in the
  others": the corrected peaks are 4.60 and 4.07 in males against 2.93 and 2.50
  in females, so two of the four are below 3.
- **Marginal-analysis sentence rests on the partition.** Averaging the four
  condition curves recovers 25% by construction, which is arithmetic rather than
  a result; the defensible statement is that a main-effect analysis captures the
  shared component and discards the 7–10% that is sex- or diet-specific.

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
