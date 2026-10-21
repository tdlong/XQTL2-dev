# Supplementary Result Y — predicting an experiment we have not run

## Methods

**Marker-based heritability.** The heritability of each trait that founder haplotypes capture genome-wide is the sum of $\hat{h}^2$ over the 268 1 cM tiles, corrected for one locus changing founder frequencies at every tile linked to it and so being counted more than once. The number of tiles a locus is counted in is the integral of its footprint, $\hat{h}^2$ against distance from the peak, divided by the peak height. We measured it at each peak of Table S1 in the treatment where that peak is strongest, on the side free of another peak within 6 cM and above the floor of $\hat{h}^2$ at non-significant tiles. Over the seventeen peaks that can be measured it ranges from 1.0 to 8.3 tiles with quartiles 2.6 and 4.1, and pooled as total integral over total height it is 3.8; at 3L 9.30 alone it is 4.1. Marker-based heritability is therefore

$$
\frac{1}{3.8}\sum_{\mathrm{tiles}} \hat{h}^2
$$

and is 11.4, 17.1, 13.4 and 16.4% for SY10 females, SY10 males, SY20 females and SY20 males. The uncertainty is in the divisor rather than in sampling: the jackknife standard error over replicates is 0.7–1.2, while taking the footprint from each peak in turn gives interquartile ranges of 10.6–16.7, 15.9–25.2, 12.5–19.7 and 15.3–24.2% (grey band, Fig. Yb). The estimate is the counterpart of SNP-based heritability at the resolution of founder haplotype blocks and, like LD score regression, applies one correction along the genome. The eighteen peaks of Table S1 sum to 72–80% of it, an in-sample upper bound on what peaks explain.

**Out-of-sample prediction.** We asked if changes in founder haplotype frequencies at leading factors can predict an unseen experiment in a manner analogous to polygenic risk scores. Each of the ten independent replicates was held out in turn and the other nine were used to build the predictor; we write ⓡ for the held-out replicate. From the other nine we estimated founder frequency changes per unit of selection, $e_{fr} = \Delta p_{fr}/i_r$ as in the heritability estimate above, and their mean $\bar{e}_f$ (the bar denotes the average over the nine non-held-out replicates in the current context); ranked tiles by heritability explained; and kept the strongest predictors separated by at least 5 cM. The predicted change in founder frequency in the held-out replicate at a given location is its known selection intensity times that mean,

$$
\widehat{\Delta p}_{fⓡ} = i_{ⓡ}\, \bar{e}_f .
$$

Heritability explained in the held-out replicate at a tile is by extension

$$
\hat{h}^2_{ⓡ} = 100\,k \sum_f \frac{\bar{e}_f\, e_{fⓡ}}{\hat{p}_f} .
$$

The $s_f^2/R$ term of the earlier heritability estimator does not appear, as $\bar{e}_f$ and $e_{fⓡ}$ are independent estimates of the same quantity, so their product has expectation $\bar{e}_f^{\,2}$. The heritability a predictor built on $N$ tiles explains in the held-out replicate (Fig. Yb) is the sum of $\hat{h}^2_{ⓡ}$ over its tiles multiplied by the ratio of that sum to the same sum with $e_{fⓡ}$ replaced by $\bar{e}_f$: the fraction of the training signal that reproduces, which is the reliability of the predictor's weights. Tiles 5 cM apart are still correlated (0.25–0.34, measured as the correlation of founder frequency change between tiles at that distance), so both sums are quadratic forms with the inverse of the $N \times N$ correlation matrix among the predictor's tiles; this changes them by at most 5%. Each held-out replicate is drawn on its own in Fig. Yb, and the ten are combined as the ratio of their means.

## Results and discussion — one paragraph

We can predict heritable variation in longevity in held-out replicates. Holding out each of the ten replicates in turn, founder-frequency changes at the *N* best markers in the other nine were used to predict the held-out replicate. For the ten best markers the prediction has *R*² = 0.60–0.69 across treatments (Fig. Ya). Marker-based heritability for the whole genome is 11–17% across treatments (Methods), of which tiles that never reach an LWP of 7.5 carry 40–47%. A predictor built on the ten best markers explains 6.4–9.5% of phenotypic variance and one built on thirty explains 8.2–11.6% (Fig. Yb), or 66–72% of the marker-based heritability. That is, thirty haplotype blocks chosen in nine experiments account for most of the heritability detectable in a tenth. Despite variation in longevity being polygenic and spread throughout the genome, a few dozen segments explain a large proportion of heritable variation.
