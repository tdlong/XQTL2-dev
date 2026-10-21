# Supplementary Result Y — predicting an experiment we have not run

## Methods

**Marker-based heritability.** The heritability of each trait that founder haplotypes capture genome-wide is the sum of $\hat{h}^2$ over the 268 1 cM tiles, corrected by a multiplicative constant to account for adjacent tiles being correlated. We estimate the constant from the peaks of Table S1. For each peak, in the treatment where it is strongest, we plotted $\hat{h}^2$ against distance from the peak, out to 6 cM on whichever side has no other peak, subtracted the floor of $\hat{h}^2$ at non-significant tiles, and divided the area under the curve by the height of the peak. Over the seventeen peaks that can be measured the signal spans from 1.0 to 8.3 tiles with a mean correction of 3.8 and quartiles of 2.6 and 4.1. Marker-based heritability is therefore

$$
\frac{1}{3.8}\sum_{\mathrm{tiles}} \hat{h}^2 .
$$

Resulting point estimates of marker-based heritability are 11.4, 17.1, 13.4 and 16.4% for SY10 females, SY10 males, SY20 females and SY20 males respectively. Dividing by the quartiles instead of 3.8 gives 10.6–16.7, 15.9–25.2, 12.5–19.7 and 15.3–24.2%. The jackknife standard error over replicates is 0.7–1.2, so the correction, not sampling, sets the uncertainty.

**Out-of-sample prediction.** We asked if changes in founder haplotype frequencies at leading factors can predict an unseen experiment in a manner analogous to polygenic risk scores. Each of the ten independent replicates was held out in turn and the other nine were used to build the predictor; we write ⓡ for the held-out replicate. From the other nine we estimated founder frequency changes per unit of selection, $e_{fr} = \Delta p_{fr}/i_r$ as in the heritability estimate above, and their mean $\bar{e}_f$ (the bar denotes the average over the nine non-held-out replicates in the current context). We ranked steps by the heritability explained in those nine and kept the best $N$, each at least 5 cM from any better one. The predicted change in founder frequency in the held-out replicate at a chosen step is its known selection intensity times that mean,

$$
\widehat{\Delta p}_{fⓡ} = i_{ⓡ}\, \bar{e}_f ,
$$

which is plotted against the observed change in Figure Ya. Heritability explained in the held-out replicate at a step is by extension

$$
\hat{h}^2_{ⓡ} = 100\,k \sum_f \frac{\bar{e}_f\, e_{fⓡ}}{\hat{p}_f} .
$$

The $s_f^2/R$ term of the earlier heritability estimator does not appear, as $\bar{e}_f$ and $e_{fⓡ}$ are independent estimates of the same quantity, so their product has expectation $\bar{e}_f^{\,2}$. Summing $\hat{h}^2_{ⓡ}$ over the $N$ steps of a predictor gives the covariance between what it predicts and what the held-out replicate does, and the same sum with $e_{fⓡ}$ replaced by $\bar{e}_f$ gives the predictor's own variance. Call these $A$ and $B$. The heritability the predictor explains in the held-out replicate is then $A^2/B$, the sum of squares accounted for by regressing observed changes on predicted ones, and $A/B$ is the reliability of weights taken from nine replicates: 0.83–0.88 for ten steps and 0.73–0.75 for thirty. Both sums are quadratic forms with the inverse of the $N \times N$ correlation matrix of the chosen steps, which corrects for their signals being correlated as the multiplicative constant does genome-wide; at a separation of 5 cM the correlation is 0.25–0.34 and the correction changes $A^2/B$ by at most 5%. Figure Yb plots $A^2/B$ against $N$, for each held-out replicate on its own and for the ten together, averaging $A$ and $B$ over replicates before forming the ratio.

## Results and discussion — one paragraph

We can predict heritable variation in longevity in held-out replicates. Holding out each of the ten replicates in turn, founder-frequency changes at the *N* best markers in the other nine were used to predict the held-out replicate. For the ten best markers the prediction has *R*² = 0.60–0.69 across treatments (Fig. Ya). Marker-based heritability for the whole genome is 11–17% across treatments (Methods), of which tiles that never reach an LWP of 7.5 carry 40–47%. A predictor built on the ten best markers explains 6.4–9.5% of phenotypic variance and one built on thirty explains 8.2–11.6% (Fig. Yb), or 66–72% of the marker-based heritability. That is, thirty haplotype blocks chosen in nine experiments account for most of the heritability detectable in a tenth. Despite variation in longevity being polygenic and spread throughout the genome, a few dozen segments explain a large proportion of heritable variation.
