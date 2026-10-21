# Supplementary Result Y — predicting an experiment we have not run

## Methods

**Out-of-sample prediction.** We asked if changes in founder haplotype frequencies at leading factors can predict an unseen experiment in a manner analogous to polygenic risk scores. Each of the ten independent replicates was held out in turn and the other nine were used to build the predictor. From those nine we estimated founder frequency changes per unit of selection, $e_{fr} = \Delta p_{fr}/i_r$ as in the heritability estimate above, ranked tiles by heritability explained, and kept the strongest predictors separated by at least 5 cM; these are the markers of Fig. Y. The predicted change in founder frequency at a chosen tile is the held-out replicate's known selection intensity, $i_r$, times the nine-replicate mean of $e_{fr}$.

The same split gives heritability explained. For held-out replicate $r$, let $\bar{e}_f$ be the mean of $e_{fr}$ over the nine training replicates and $e_{fr}$ the held-out replicate's own estimate; heritability explained at a tile is

$$
\hat{h}^2_r = 100\,k \sum_f \frac{\bar{e}_f\, e_{fr}}{\hat{p}_f} .
$$

Squaring one noisy estimate inflates it by its sampling variance, which the $s_f^2/R$ term of the full estimator removes; $\bar{e}_f$ and $e_{fr}$ are independent estimates of the same quantity, so their product has expectation $\bar{e}_f^{\,2}$ with no such term, and as tiles were chosen on the training replicates alone the estimate carries no selection bias. Fig. Yb reports, for the $N$ highest-ranked tiles, the heritability the predictor explains: the sum of $\hat{h}^2_r$ over those tiles multiplied by its ratio to the training value $100k\sum_f \bar{e}_f^{\,2}/\hat{p}_f$, the fraction of the training heritability recovered in the held-out replicate.

Tiles in linkage share signal, so summing over tiles counts a locus once for every tile it reaches. For the $N$ tiles of a predictor we combine them through $C$, the $N \times N$ matrix of correlations among them, estimated as the correlation of $\bar{e}$ at one tile with $e_r$ at another over the ten folds: writing $\tilde{e}_{fr} = e_{fr}/\sqrt{\hat{p}_f}$ as a vector over the $N$ tiles, the sum above becomes $100k\sum_f \tilde{\bar{e}}_f^{\top} C^{-1} \tilde{e}_{fr}$; at 5 cM spacing $C$ is near the identity and the correction is about 1%. Marker-based heritability of the whole genome (Fig. Yb, dashed) is the sum of $\hat{h}^2_r$ over all 273 tiles divided by the mean number of tiles a locus is counted in, $\operatorname{tr}(C^2)/273 = 4.4$–$5.1$, with $C$ built from the same correlation at increasing tile separation (0.88 at 1 cM, 0.25 at 5 cM, zero by 15 cM). It is the counterpart of SNP-based heritability and, like LD score regression, assumes one correction along the genome. Jackknife over the ten folds gives standard errors of 2.6–4.1 percentage points.

## Results and discussion — one paragraph

We can predict heritable variation in longevity in held-out replicates. Holding out each of the ten replicates in turn, founder-frequency changes at the *N* best markers in the other nine were used to predict the held-out replicate. For the ten best markers the prediction has *R*² = 0.60–0.69 across treatments (Fig. Ya). Marker-based heritability for the whole genome is 10–14% across treatments (SE 2.6–4.1; Methods). A predictor built on the ten best markers explains 6.4–9.5% of phenotypic variance and one built on thirty explains 8.2–11.6% (Fig. Yb), 82–91% of the marker-based heritability. Thirty haplotype blocks chosen in nine experiments account for most of the heritability detectable in a tenth: variation in longevity is spread across the genome, but across a few dozen segments whose effects reproduce, not a continuum of effects too small to find.
