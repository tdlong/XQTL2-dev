# Supplementary Result Y — predicting an experiment we have not run

## Methods

**Marker-based heritability.** The heritability of each trait that founder haplotypes capture genome-wide is the sum of $\hat{h}^2$ over the 268 1 cM tiles, corrected by a multiplicative constant to account for adjacent tiles being correlated. We estimate the constant from the peaks of Table S1. For each peak, in the treatment where it is strongest, we plotted $\hat{h}^2$ against distance from the peak, out to 6 cM on whichever side has no other peak, subtracted the floor of $\hat{h}^2$ at non-significant tiles, and divided the area under the curve by the height of the peak. Over the seventeen peaks that can be measured the signal spans from 1.0 to 8.3 tiles with a mean correction of 3.8 and quartiles of 2.6 and 4.1. Marker-based heritability is therefore

$$
\frac{1}{3.8}\sum_{\mathrm{tiles}} \hat{h}^2 .
$$

Resulting point estimates of marker-based heritability are 11.4, 17.1, 13.4 and 16.4% for SY10 females, SY10 males, SY20 females and SY20 males respectively. Dividing by the quartiles instead of 3.8 gives 10.6–16.7, 15.9–25.2, 12.5–19.7 and 15.3–24.2%.

**Out-of-sample prediction.** A polygenic risk score predicts a phenotype from a genotype it played no part in measuring. The counterpart here is to take the control pool of a held-out replicate as given — its founder frequencies are observed without any selection — and predict the founder frequencies of its *selected* pool, scoring the prediction against what that pool actually shows. Each of the ten replicates was held out in turn; we write ⓡ for the held-out one. From the other nine we estimated the response per unit of selection, $e_{fr} = \Delta p_{fr}/i_r$ as in the heritability estimate above, with mean $\bar{e}_f$; ranked positions by the heritability explained in those nine; and kept the best $N$ as markers, each at least 5 cM from any better one. Since the held-out replicate's selected proportion is known, its predicted selected-pool frequency at a marker is

$$
\hat{z}_{fⓡ} = \hat{c}_{fⓡ} + i_{ⓡ}\, \bar{e}_f ,
$$

with $\hat{c}_{fⓡ}$ the observed control frequency. Figure Ya plots the predicted change, $\hat{z}_{fⓡ} - \hat{c}_{fⓡ}$, against the observed change rather than the frequencies themselves: the control frequency enters both and would inflate the agreement.

Heritability explained in the held-out replicate at a marker is by extension

$$
\hat{h}^2_{ⓡ} = 100\,k \sum_f \frac{\bar{e}_f\, e_{fⓡ}}{\hat{p}_f} .
$$

The $s_f^2/R$ term of the earlier heritability estimator does not appear, as $\bar{e}_f$ and $e_{fⓡ}$ are independent estimates of the same quantity, so their product has expectation $\bar{e}_f^{\,2}$. Markers in linkage tag the same region and their heritabilities cannot simply be added, the double counting the constant above corrects genome-wide; among the markers of a predictor it is removed by inverting their measured correlation, which at a separation of 5 cM changes the sum by at most 5%. Two quantities follow, and they answer different questions. The first is what the predictor explains: the sum of $\hat{h}^2_{ⓡ}$ over the $N$ markers, squared and divided by the same sum with $e_{fⓡ}$ replaced by $\bar{e}_f$, which is the squared covariance between predicted and observed changes over the variance of the predicted ones — a predictor's share of variance, as in any polygenic score. The second is how much heritability the $N$ regions hold, whoever estimates their effects. Ranking and effect estimation must then be separated, or a region is credited using the data that selected it: five replicates rank the positions, four others estimate $\bar{e}_f$, and the held-out replicate scores the result, averaged over five random splits of the nine. Without that separation the regions appear to hold about 10% more than they do.

## Results and discussion — one paragraph

We can predict heritable variation in longevity in held-out replicates. Holding out each of the ten replicates in turn, its control pool and the frequency changes at the *N* best markers in the other nine predict the founder frequencies of its selected pool; for the ten best markers, predicted and observed changes agree with *R*² = 0.60–0.69 across treatments (Fig. Ya). A predictor built on thirty markers explains 8.3–12.7% of phenotypic variance in the held-out replicate, against a marker-based heritability of 11.4–17.1% for the whole genome (Fig. Yb; Methods). The thirty regions themselves hold 9.6–13.9%. A predictor therefore captures 84–91% of what its regions hold, while those regions hold 78–84% of the heritability in the genome: the limit is the number of regions, not the precision with which their effects are estimated. Both curves are still rising at forty markers, where 5 cM clumping exhausts the genome, and the gap between them stays small. Despite variation in longevity being polygenic and spread throughout the genome, a few dozen haplotype blocks account for most of the heritable variation we can detect.

## Figure Y legend

**Figure Y: Out of sample predictions.** A, Founder frequency change at the top 10 markers in nine replicates predicts the change in a tenth, held-out replicate. Colors are founder haplotypes as in Figure 2. The dashed line is the line of identity. B, Heritability as a function of the number of markers used, in percent of phenotypic variance. Red, the heritability the *N* regions hold, with the regions chosen in five replicates and their effects estimated in four others, so that no region is credited using the data that selected it. Blue, the heritability an *N*-marker predictor built on nine replicates explains in the held-out replicate; light lines are the ten replicates individually. The dashed line is the marker-based heritability of the whole genome and the grey band its interquartile range as the footprint width is taken from each peak of Table S1 in turn.
