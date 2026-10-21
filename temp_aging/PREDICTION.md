# Supplementary Result Y — predicting an experiment we have not run

## Methods

**Marker-based heritability.** The heritability of each trait that founder haplotypes capture genome-wide is the sum of $\hat{h}^2$ over the 268 1 cM tiles, corrected for one locus changing founder frequencies at every tile linked to it and so being counted more than once. The number of tiles a locus is counted in is the integral of its footprint, $\hat{h}^2$ against distance, divided by its height. Over the peaks of Table S1, each in its strongest treatment and on the side free of another peak within 6 cM, $\hat{h}^2$ falls to a tenth of the peak by 2–6 cM and the footprint integral is 3.8 tiles pooled over peaks (median 3.4; 4.1 at 3L 9.30). Marker-based heritability is therefore

$$
\frac{1}{3.8}\sum_{\mathrm{tiles}} \hat{h}^2 ,
$$

and is 11.4, 17.1, 13.4 and 16.4% for SY10 females, SY10 males, SY20 females and SY20 males, with jackknife standard errors over replicates of 0.7–1.2. It is the counterpart of SNP-based heritability at the resolution of founder haplotype blocks, and like LD score regression assumes one correction factor along the genome. The eighteen peaks of Table S1 sum to 72–80% of it, an in-sample upper bound on what peaks explain.

**Out-of-sample prediction.** We asked if changes in founder haplotype frequencies at leading factors can predict an unseen experiment in a manner analogous to polygenic risk scores. Each of the ten independent replicates was held out in turn and the other nine were used to build the predictor; we write ⓡ for the held-out replicate. From the other nine we estimated founder frequency changes per unit of selection, $e_{fr} = \Delta p_{fr}/i_r$ as in the heritability estimate above, and their mean $\bar{e}_f$ (the bar denotes the average over the nine non-held-out replicates in the current context); ranked tiles by heritability explained; and kept the strongest predictors separated by at least 5 cM. The predicted change in founder frequency in the held-out replicate at a given location is its known selection intensity times that mean,

$$
\widehat{\Delta p}_{fⓡ} = i_{ⓡ}\, \bar{e}_f .
$$

Heritability explained in the held-out replicate at a tile is by extension

$$
\hat{h}^2_{ⓡ} = 100\,k \sum_f \frac{\bar{e}_f\, e_{fⓡ}}{\hat{p}_f} .
$$

The $s_f^2/R$ term of the earlier heritability estimator does not appear, as $\bar{e}_f$ and $e_{fⓡ}$ are independent estimates of the same quantity, so their product has expectation $\bar{e}_f^{\,2}$. The tiles of a predictor are at least 5 cM apart, where the correlation between tiles is 0.25–0.34, and their heritabilities are added without correction; dividing out the $N \times N$ correlation among them changes the sums by at most 5%. The heritability a predictor built on $N$ tiles explains in the held-out replicate (Fig. Yb) is the sum of $\hat{h}^2_{ⓡ}$ over its tiles multiplied by the ratio of that sum to the same sum with $e_{fⓡ}$ replaced by $\bar{e}_f$: the fraction of the training signal that reproduces, which is the reliability of the predictor's weights. The dashed line in Fig. Yb is the marker-based heritability defined above.

## Results and discussion — one paragraph

We can predict heritable variation in longevity in held-out replicates. Holding out each of the ten replicates in turn, founder-frequency changes at the *N* best markers in the other nine were used to predict the held-out replicate. For the ten best markers the prediction has *R*² = 0.60–0.69 across treatments (Fig. Ya). Marker-based heritability for the whole genome is 11–17% across treatments (SE 0.7–1.2; Methods), of which tiles that never reach an LWP of 7.5 carry 40–47%. A predictor built on the ten best markers explains 6.4–9.5% of phenotypic variance and one built on thirty explains 8.2–11.6% (Fig. Yb), or 66–72% of the marker-based heritability. That is, thirty haplotype blocks chosen in nine experiments account for most of the heritability detectable in a tenth. Despite variation in longevity being polygenic and spread throughout the genome, a few dozen segments explain a large proportion of heritable variation.
