# Supplementary Result Y — predicting an experiment we have not run

## Methods

**Out-of-sample prediction.** We asked if changes in founder haplotype frequencies at leading factors can predict an unseen experiment, by analogy with polygenic risk scores. Each of the ten independent replicate cages was held out in turn and the other nine were used to build the predictor. From those nine we estimated founder frequency changes (the mean of the $e_{fr}$ defined above), ranked tiles by heritability explained, and kept the strongest predictors separated by at least 5 cM. The predicted change in founder frequency at chosen positions is simply the known selection intensity, $i$, times the nine-cage estimate. We similarly estimate heritability explained by the estimator.

At each chosen tile heritability explained is estimated from the two halves of the split, the nine-cage mean $\bar{e}_{f,\mathrm{nine}}$ and the held-out cage's $e_{f,\mathrm{tenth}}$, as

$$
\hat{h}^2 = 100\,k \sum_f \frac{\bar{e}_{f,\mathrm{nine}}\; e_{f,\mathrm{tenth}}}{\hat{p}_f}.
$$

The two factors have independent sampling errors, so no $s_f^2/R$ term is required, and because tiles were chosen on the nine cages alone the estimate is not inflated by that selection. The heritability a score on $N$ tiles explains (Fig. Yb) is this quantity scaled by the fraction of the nine-cage heritability, $100k\sum_f \bar{e}_{f,\mathrm{nine}}^{\,2}/\hat{p}_f$, that is recovered in the tenth cage.

Summing over tiles counts a locus once for every tile it is linked to. For the $N$ tiles of a score we correct for this with the $N \times N$ matrix of correlations among them, estimated from the same two halves (the nine-cage response at one tile against the tenth-cage response at another); at 5 cM spacing the matrix is close to diagonal and the correction is about 1%. The same quantity over all 273 tiles, the analogue of SNP-based heritability, we report only as an order of magnitude: adjacent 1 cM tiles are 89% correlated, the genome comprises about 55 independent segments, and across defensible choices of correction the total ranges from 6 to 32%.

## Results and discussion — one paragraph

Haplotype frequency change in a population we have never seen is predictable from
the other nine. At the ten highest-ranked markers, predicted and observed changes
correlate at $r = 0.77$–$0.83$ across the four treatments, with regression slopes
of 0.87–0.91 against a line of identity that was not fitted (Fig. Y a); the
prediction is of eighty signed numbers and is as well calibrated as it is
correlated, founder B3 predicted to fall 11.4% at chr3L 9.31 Mb in one held-out
cage and falling 12.0%. A score built on those ten markers explains 6.4–9.5% of
phenotypic variance, and one built on thirty explains 8.2–11.6% (Fig. Y b),
against a marker-based heritability for the whole genome of order 10%. This is
the validation a selection experiment can offer in place of functional work: the
mapped signal predicts quantitatively, genome-wide, in material that contributed
nothing to the mapping. It does not place causal sites within the mapped
intervals — a linked haplotype predicts a held-out population as well as a causal
one — nor test whether effects combine additively across loci, which pooled
sequencing cannot see. What limits the genome-wide total is the same feature that
makes the per-marker prediction precise: eight founders at appreciable frequency
concentrate the segregating variation into roughly 55 independent segments, each
assayed in thousands of flies. A human association study has a thousandfold more
resolution and no comparable precision at a locus, and the two cannot be had at
once.
