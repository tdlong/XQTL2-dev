# Supplementary Result Y — predicting an experiment we have not run

## Methods

**Leave one cage out.** The ten replicate cages are the unit. Holding out cage $r$, the other nine give a standardised response per founder at every step, $e_f = \overline{\Delta p_f / i}$, where $i$ is the selection intensity the replicate's selected proportion implies. The held-out cage's own proportion is known, so its predicted founder frequency change is $i_r e_f$ with nothing fitted to it. Positions are ranked by the heritability of the nine training cages and greedily clumped 5 cM apart, so the held-out cage neither sets the weights nor chooses where they are read. Prediction is scored over all ten folds and all four treatments.

**An unbiased heritability.** The estimator of the main text squares a replicate's response and subtracts $s_f^2/R$ to remove the inflation squaring introduces. Training and held-out cages are different flies, so their sampling noise is independent and the cross-product needs no such correction,

$$
\hat{h}^2 = 100\,k \sum_f \frac{e_f^{\,\mathrm{train}}\, e_f^{\,\mathrm{test}}}{\hat{p}_f}
$$

since $E[e^{\mathrm{train}} e^{\mathrm{test}}] = e^2$. Because positions are chosen from the training cages alone it also carries no winner's curse: a position selected on noise contributes zero in expectation rather than a positive bias.

**Correlation between markers.** Heritabilities at nearby positions cannot be added, because one locus raises every position it is linked to. We measure that correlation as the cross-fold product at lag $d$ — the training response at one tile against the held-out response $d$ cM away — which carries no noise contribution at any lag, including zero, where the within-sample autocorrelation is inflated by the $\pm100$ kb smoothing. The product between arms, which cannot reflect linkage, is $-1.9\%$ of the lag-zero value, so there is no genome-wide common term to remove. The correlation falls from 0.88 at 1 cM to 0.25 at 5 cM and reaches zero by 15 cM.

**Marker-based heritability.** Writing $u_f = e_f/\sqrt{\hat{p}_f}$ and $R$ for the correlation between markers, the heritability a set of markers explains is the quadratic form $100\,k \sum_f u_f^{\mathrm{train}\top} R^{-1} u_f^{\mathrm{test}}$, not the sum of the separate values. This is what GREML estimates through the relationship matrix, so it is the synthetic-population counterpart of SNP-based heritability, and is reported as such rather than as the trait's heritability.

Among the markers of a score $R$ is small and well conditioned, 1.4 to 3.1 at these spacings, and the correction is slight: naive and corrected sums agree to about 1%, so markers 5 cM apart are already close to independent. Genome-wide the same form over all 273 one-cM tiles cannot be used, because adjacent tiles correlate at 0.88 and the matrix is near-singular. Thinning the tiles to the correlation length restores conditioning, and over 4–5 cM spacing — where the matrix is well conditioned and few markers have yet been discarded — it agrees with the scalar approximation $\operatorname{tr}(R^2)/n$, which requires no inverse and is stable to every choice entering it. We report the scalar estimate.

**Why this may approach the trait's heritability.** In human studies SNP-based heritability falls well short of the trait's, because common markers tag causal variation imperfectly and rare variants are missed. Neither limitation operates here. Eight founder haplotypes exhaust the variation segregating in the population, their frequencies are estimated directly rather than tagged, and any variant differing among founders is therefore captured exactly. What remains outside is variation not segregating among the eight founders, and regions where founders are locally indistinguishable — small terms, where in human studies the gap is a factor of two. Marker-based and narrow-sense heritability should therefore be close in this population, which is not an assumption a human study can make.

## Results

Nine cages predict the tenth. Over the ten highest-ranked positions, predicted and observed founder frequency changes in a cage that contributed to neither the weights nor the choice of position correlate at $r = 0.77$ to $0.83$, with a slope of 0.87 to 0.91 on a line of identity that was not fitted (Supplementary Figure Y, panel A). The prediction is of eighty signed numbers and is as well calibrated as it is correlated: founder B3 was predicted to fall 11.4% at chr3L 9.31 Mb in a held-out cage and fell 12.0%.

Marker-based heritability of the whole genome is 10.0, 14.0, 10.7 and 12.0% of phenotypic variance in the four treatments, with jackknife standard errors of 2.6 to 4.1, against naive sums over tiles of 44 to 64% that the correlation between tiles accounts for.

Against that, a score of a few markers explains most of it.

**(a) Heritability carried by the top $N$ markers** (% of phenotypic variance)

| $N$ | SY10 F | SY10 M | SY20 F | SY20 M |
|---:|---:|---:|---:|---:|
| 5 | 5.7 | 8.1 | 5.4 | 8.3 |
| 10 | 7.7 | 11.2 | 8.4 | 10.7 |
| 15 | 9.0 | 12.7 | 9.9 | 11.8 |
| 20 | 9.8 | 13.9 | 11.1 | 13.0 |
| 30 | 11.2 | 15.5 | 12.1 | 14.6 |

**(b) Fraction of it the predictor captures**

| $N$ | SY10 F | SY10 M | SY20 F | SY20 M |
|---:|---:|---:|---:|---:|
| 5 | 0.88 | 0.86 | 0.89 | 0.89 |
| 10 | 0.83 | 0.85 | 0.88 | 0.86 |
| 15 | 0.79 | 0.81 | 0.83 | 0.83 |
| 20 | 0.76 | 0.80 | 0.81 | 0.81 |
| 30 | 0.73 | 0.75 | 0.73 | 0.75 |

**(c) Heritability the score explains, (a) × (b)** (% of phenotypic variance)

| $N$ | SY10 F | SY10 M | SY20 F | SY20 M |
|---:|---:|---:|---:|---:|
| 5 | 5.0 | 7.0 | 4.8 | 7.3 |
| 10 | 6.4 | 9.5 | 7.4 | 9.2 |
| 15 | 7.1 | 10.2 | 8.2 | 9.8 |
| 20 | 7.4 | 11.1 | 8.9 | 10.5 |
| 30 | 8.2 | 11.6 | 8.9 | 10.9 |

Ten markers explain 6.4 to 9.5% of phenotypic variance, 64 to 77% of the marker-based heritability; thirty explain 8.2 to 11.6%, or 82 to 91% (Supplementary Figure Y, panel B). The fraction captured falls as markers are added, from 0.88 at five to 0.74 at thirty, because later markers carry weaker signal and noisier weights, while the heritability at the markers continues to rise, so their product still climbs.

This substitutes for functional validation in the one respect functional validation would also address — that the mapped signal is real rather than a threshold artefact — and does so genome-wide rather than one locus at a time. It does not establish that the mapped intervals contain the causal sites: a linked haplotype predicts a held-out cage exactly as well as a causal one. Nor does it test whether effects combine additively across loci, which pooled sequencing cannot see, a pool reporting only the mean of its individuals' genotypes.
