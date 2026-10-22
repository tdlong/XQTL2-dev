# Supplementary Result Y — predicting an experiment we have not run

## Methods

**Marker-based heritability.** The heritability of each trait that founder haplotypes capture genome-wide is the sum of $\hat{h}^2$ over the 268 1 cM tiles, corrected by a multiplicative constant to account for adjacent tiles being correlated. We estimate the constant from the peaks of Table S1. For each peak, in the treatment where it is strongest, we plotted $\hat{h}^2$ against distance from the peak, out to 6 cM on whichever side has no other peak, subtracted the floor of $\hat{h}^2$ at non-significant tiles, and divided the area under the curve by the height of the peak. Over the seventeen peaks that can be measured the signal spans from 1.0 to 8.3 tiles with a mean correction of 3.8 and quartiles of 2.6 and 4.1. Marker-based heritability is therefore

$$
\frac{1}{3.8}\sum_{\mathrm{tiles}} \hat{h}^2 .
$$

Resulting point estimates of marker-based heritability are 11.4, 17.1, 13.4 and 16.4% for SY10 females, SY10 males, SY20 females and SY20 males respectively. Dividing by the quartiles instead of 3.8 gives 10.6–16.7, 15.9–25.2, 12.5–19.7 and 15.3–24.2%.

**Out-of-sample prediction.** A polygenic risk score predicts a phenotype from an out of sample individual. The counterpart here is to take the control pool genotypes of a held-out replicate as given and predict the founder frequencies of its *selected* pool, then score those predictions against what that pool shows. Each of the ten replicates was held out in turn and we write ⓡ to indicate the held-out sample below. From the other nine we estimate the response per unit of selection, $e_{fr} = \Delta p_{fr}/i_r$ as in the heritability estimate above, with mean $\bar{e}_f$; rank positions by the heritability explained in those nine; and keep the best $N$ as markers, each at least 5 cM from any better one. Since the held-out replicate's selected proportion is known, its predicted selected-pool frequency at a marker is

$$
\hat{z}_{fⓡ} = \hat{c}_{fⓡ} + i_{ⓡ}\, \bar{e}_f ,
$$

with $\hat{c}_{fⓡ}$ the observed control frequency. Figure Ya plots the predicted change, $\hat{z}_{fⓡ} - \hat{c}_{fⓡ}$, against the observed changes.

Heritability explained in the held-out replicate at a marker is by extension

$$
\hat{h}^2_{ⓡ} = 100\,k \sum_f \frac{\bar{e}_f\, e_{fⓡ}}{\hat{p}_f} .
$$

The $s_f^2/R$ term of the earlier heritability estimator does not appear, as $\bar{e}_f$ and $e_{fⓡ}$ are independent estimates of the same quantity, so their product has expectation $\bar{e}_f^{\,2}$. Markers in linkage tag the same region and their heritabilities cannot simply be added. We therefore summed $\hat{h}^2_{ⓡ}$ over the $N$ markers of a predictor using the inverse of their correlation matrix, estimated from the correlation in founder frequency change between markers as a function of distance. Because markers are at least 5 cM apart, this changes the sum over simple addition by at most 5%.

We used this estimate of heritability explained in two ways. First, to estimate the heritability explained by the best $N$ markers themselves (red line, Fig. Yb). If the same replicates both choose the markers and estimate their effects, the estimate is inflated by the winner's curse. We therefore split the nine non-held-out replicates: five were used to rank positions and choose the $N$ markers, the other four to estimate $\bar{e}_f$, and the held-out replicate to estimate the heritability explained, averaging over five random splits. Omitting the split inflates the estimate by 12–17%. This term can be thought of as a perfect best N-marker estimator. Second, to estimate the heritability explained by an $N$-marker predictor whose effects come from all nine replicates (blue line). This is the squared covariance between predicted and observed frequency changes divided by the variance of the predicted changes. That is, the sum of $\hat{h}^2_{ⓡ}$ over the $N$ markers, squared, divided by the same sum with $e_{fⓡ}$ replaced by $\bar{e}_f$. It is smaller than the heritability explained by the markers themselves because $\bar{e}_f$ is an estimate, and can be thought of as how well we predict heritability in a heldout sample where founder allele frequencies are unknown in the selected pool.

## Results and discussion — one paragraph

We can predict heritable variation in longevity in held-out replicates. Holding out each of the ten replicates in turn, its control pool and the frequency changes at the *N* best markers in the other nine predict the founder frequencies of its selected pool. For the ten best markers, predicted and observed founder frequency changes agree with *R*² = 0.60–0.69 across treatments (Fig. Ya). Predicting founder frequencies in held-out selected pools further allows us to estimate the proportion of heritable variation in longevity explained. We first asked how much heritable variation the best *N* markers explain. To do this we chose the markers in five replicates, estimated their effects in four others, and estimated the heritability they explain in the tenth. This sampling scheme accounts for the winner's curse in choosing the best *N* markers. The best 20–30 markers explain 71–84% of marker-based heritability (red line, Fig. Yb), with the remainder presumably in smaller-effect polygenic regions. We then asked what fraction of this maximum is captured by a predictor whose effects are estimated from all nine replicates. At thirty markers the predictor explains 8.3–12.7% of phenotypic variance, which is 84–91% of what those markers explain. The predictor loses efficiency as markers are added, apparent as the blue lines beginning to plateau at around 20–30 markers. This is presumably because effects at more polygenic loci are estimated less precisely, so adding them adds little additional predictive value. Caveats aside, despite variation in longevity being polygenic and spread throughout the genome, a few dozen haplotype blocks predictably account for a large fraction of the heritable variation we can detect.

## Figure Y legend

**Figure Y: Out of sample predictions.** A, Founder frequency change at the top 10 markers in nine replicates predicts the change in a tenth, held-out replicate. Colors are founder haplotypes as in Figure 2. The dashed line is the line of identity. B, Heritability as a function of the number of markers used, as a percent of phenotypic variation. The dashed line is the marker-based heritability of the whole genome and the grey band its interquartile range as the footprint width is taken from each peak of Table S1 in turn. The red line is the heritability explained by the best *N* markers, with the markers chosen in five replicates, their effects estimated in four others, and the heritability they explain estimated in the held-out replicate. The blue line is the heritability an *N*-marker predictor built on nine replicates explains in the held-out replicate, with the lighter blue lines each individual replicate.
