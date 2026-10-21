# Supplementary Result Y — predicting an experiment we have not run

## Methods

**Out-of-sample prediction.** We asked if changes in founder haplotype frequencies at leading factors can predict an unseen experiment, by analogy with polygenic risk scores. Each of the ten independent replicate cages was held out in turn and the other nine were used to build the predictor. From those nine we estimated founder frequency changes (the mean of the $e_{fr}$ defined above), ranked tiles by heritability explained, and kept the strongest predictors separated by at least 5 cM. The predicted change in founder frequency at chosen positions is simply the known selection intensity, $i$, times the nine-cage estimate. We similarly estimate heritability explained by the estimator.

At each chosen tile the nine cages give $\bar{e}_{f,\mathrm{nine}}$, the mean of $e_{fr}$ over those cages, and the held-out cage gives $e_{f,\mathrm{tenth}}$. Heritability explained at the tile is

$$
\hat{h}^2 = 100\,k \sum_f \frac{\bar{e}_{f,\mathrm{nine}}\; e_{f,\mathrm{tenth}}}{\hat{p}_f}
$$

where $\hat{p}_f$ is the mean control frequency. The two factors carry independent sampling error, so no $s_f^2/R$ term is needed, and as tiles were chosen on the nine cages alone the estimate is free of winner's curse. The heritability a score on $N$ tiles explains (Fig. Yb) is this cross-product multiplied by its ratio to the nine-cage estimate, $100k\sum_f \bar{e}_{f,\mathrm{nine}}^{\,2}/\hat{p}_f$.

Tiles in linkage share signal, so heritability over $N$ tiles is combined through the $N \times N$ matrix $R$ of correlations among them, estimated as the correlation of $\bar{e}_{\mathrm{nine}}$ at one tile with $e_{\mathrm{tenth}}$ at another,

$$
\hat{h}^2_N = 100\,k \sum_f \mathbf{u}_{f,\mathrm{nine}}^{\top} R^{-1} \mathbf{u}_{f,\mathrm{tenth}}
$$

with $u_f = e_f/\sqrt{\hat{p}_f}$ over the $N$ tiles. At 5 cM spacing $R$ is close to the identity and the correction is about 1%. Over all 273 tiles adjacent entries of $R$ are 0.89 and the genome holds about 55 independent segments; the total, the analogue of SNP-based heritability, is then not well determined, ranging from 6 to 32% across methods, and is reported as an order of magnitude.

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
