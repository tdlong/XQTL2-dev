# Supplementary Result Y — predicting an experiment we have not run

## Methods

**Out-of-sample prediction.** We asked if changes in founder haplotype frequencies at leading factors can predict an unseen experiment, by analogy with polygenic risk scores. Each of the ten independent replicates was held out in turn and the other nine were used to build the predictor. From those nine we estimated founder frequency changes (the mean of the $e_{fr}$ defined above), ranked tiles by heritability explained, and kept the strongest predictors separated by at least 5 cM. The predicted change in founder frequency at chosen positions is simply the known selection intensity, $i$, times the estimate from the other nine replicates. We similarly estimate heritability explained by the predictor.

For held-out replicate $r$, write $\bar{e}_{f,-r}$ for the mean of $e_{fr'}$ over the nine remaining replicates. Heritability explained at a chosen tile is

$$
\hat{h}^2_r = 100\,k \sum_f \frac{\bar{e}_{f,-r}\; e_{fr}}{\hat{p}_f}
$$

where $e_{fr}$ is the held-out replicate's own response. The two factors carry independent sampling error, so the $s_f^2/R$ term is not required, and as tiles were chosen on the nine remaining replicates the estimate is not inflated by that choice. Heritability explained by a predictor built on $N$ tiles (Fig. Yb) is this cross-product multiplied by its ratio to the nine-replicate estimate $100k\sum_f \bar{e}_{f,-r}^{\,2}/\hat{p}_f$.

Tiles in linkage share signal, so heritability over the $N$ tiles is combined through $C$, the $N \times N$ matrix of correlations among them, estimated as the correlation of $\bar{e}_{-r}$ at one tile with $e_r$ at another. Writing $\tilde{e}_{fr} = e_{fr}/\sqrt{\hat{p}_f}$ as a vector over the $N$ tiles,

$$
\hat{h}^2_{N,r} = 100\,k \sum_f \tilde{e}_{f,-r}^{\top} C^{-1} \tilde{e}_{fr} .
$$

At 5 cM spacing $C$ is close to the identity and the correction is about 1%. Over all 273 tiles adjacent entries of $C$ are 0.89 and the genome comprises about 55 independent segments; the marker-based heritability of the whole genome, the counterpart of SNP-based heritability, is then not well determined, ranging from 6 to 32% across methods, and is reported as an order of magnitude only.

## Results and discussion — one paragraph

Haplotype frequency change in a population we have never seen is predictable from
the other nine. At the ten highest-ranked tiles, predicted and observed changes
correlate at $r = 0.77$–$0.83$ across the four treatments, with regression slopes
of 0.87–0.91 against a line of identity that was not fitted (Fig. Y a); the
prediction is of eighty signed numbers and is as well calibrated as it is
correlated, founder B3 predicted to fall 11.4% at chr3L 9.31 Mb in one held-out
cage and falling 12.0%. A score built on those ten tiles explains 6.4–9.5% of
phenotypic variance, and one built on thirty explains 8.2–11.6% (Fig. Y b),
against a marker-based heritability for the whole genome of order 10%. This is
the validation a selection experiment can offer in place of functional work: the
mapped signal predicts quantitatively, genome-wide, in material that contributed
nothing to the mapping. It does not place causal sites within the mapped
intervals — a linked haplotype predicts a held-out population as well as a causal
one — nor test whether effects combine additively across tiles, which pooled
sequencing cannot see. What limits the genome-wide total is the same feature that
makes the per-marker prediction precise: eight founders at appreciable frequency
concentrate the segregating variation into roughly 55 independent segments, each
assayed in thousands of flies. A human association study has a thousandfold more
resolution and no comparable precision at a locus, and the two cannot be had at
once.
