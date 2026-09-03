# Methods (incomplete)

Library preparation and sequencing: DNA was prepared on pools of several hundred flies using ?? (Supp table of N of flies).  Illumina DNAseq kit was used for library preparation and was sequenced at UCI Genomics Research and Technology Hub on an Illumina NovaSeq X Plus as 2 × 150 bp reads. Sixty experimental libraries were sequenced: ten replicates of the experiment cages × three groups (unselected control, SY10 selected, SY20 selected) × two sexes.  Reads were aligned to the Drosophila melanogaster dm6 reference genome with bwa mem (Li 2013), coordinate sorted and indexed. Experimental sample alignments were not deduplicated although reads from the eight founder alignments were.

SNP calling: SNPs were called on the five major chromosome arms for the eight founder samples and a catalog of SNPs segregating among the eight founders of the synthetic population created. Founders were called jointly with bcftools (cite): a site entered the catalog if it was biallelic, covered at 10× or more in every founder, had a minor allele frequency of 3% or less in every founder, segregating an alternate allele in at least one founder, and was at least 20 bp from the nearest founder indel. The final catalog comprises 1,207,436 SNPs from 1,887,667 initial candidate positions. For each of the 60 experimental samples reference and alternate read counts were then taken directly from the pileup at each catalog position, and these counts used to infer founder haplotype frequencies. At called SNPs autosomal depth averaged 140× per library and a range of 49× to 400×; fifteen libraries fell below 75× and twenty below 100×. In prior work we have shown 50X coverage is sufficient for high quality haplotype estimation.

Haplotype estimation: The synthetic population is derived from eight inbred founders, so a pool's allele frequency at any SNP is the sum over founders of each founder's frequency in the pool weighted by its genotype. The founder genotypes are known, so over a local window consisting of many SNPs the eight founder frequencies can be estimated via constrained least squares (lsei, R package limSolve) subject to the frequencies summing to one and each being non-negative. Those frequencies were estimated at 22,420 five kb euchromatic region *steps*, with each step given an initial *window* of ±75 kb, widened in inverse proportion to the local recombination rate so that windows span comparable genetic rather than physical distance. The fit returns a covariance matrix for the estimates which was averaged over the ten steps either side of each step and used in the Wald test (below). Founders genotypes were hierarchically clustered within each window, cut at a height of 2.5, and indistinguishable founders filled by linear interpolation between the means of the flanking resolved steps. The haplotype frequencies were smoothed with a running mean of ±100 kb (41 steps).

**Wald test.** At each step let $\hat{c}$ and $\hat{z}$ be the eight-vectors of estimated founder frequencies in the control and selected pools, and $V_c$, $V_z$ their covariance matrices. The Wald-statistic is

$$
T = \frac{1}{R^2}\,(\hat{c} - \hat{z})^{\top}(V_c + V_z)^{-}(\hat{c} - \hat{z})
$$

on seven degrees of freedom, since the eight haplotype frequencies sum to one. That constraint also makes $V_c + V_z$ singular, so $(\cdot)^{-}$ is a generalised inverse, obtained by eigendecomposition with the null direction discarded and the remaining eigenvalues floored at one per cent of their mean, so that near-degenerate directions cannot inflate $T$. Significance is reported as $-\log_{10} P$.

Each pool's covariance is $V = M + \Sigma$, where $M$ is the multinomial sampling covariance on sampled chromosomes with subscript denoting founder alleles,

$$
M_{ff} = \frac{\bar{p}_f(1-\bar{p}_f)}{n_{\mathrm{eff}}}, \qquad M_{fg} = -\frac{\bar{p}_f \bar{p}_g}{n_{\mathrm{eff}}}
$$

evaluated at the frequency $\bar{p}$ pooled across the selected and control samples, and $\Sigma$ is the smoothed haplotype reconstruction covariance returned by the constrained least squares fit above. The effective number of individuals is estimated from a pool of $n$ flies (where $k$ is one for male X chromosomes and two otherwise) as

$$
n_{\mathrm{eff}} = \frac{kn \sum_f \bar{p}_f(1-\bar{p}_f)}{\sum_f \bar{p}_f(1-\bar{p}_f) + kn \operatorname{tr}(\Sigma)}
$$

The 10 replicates were combined as an $n_{\mathrm{eff}}$-weighted mean of $\hat{c}$ and $\hat{z}$, with covariances similarly pooled. As the haplotype frequencies entering $T$ are smoothed, their variance is reduced relative to unsmoothed estimates. $R^2$, the squared correlation between raw and smoothed frequencies taken over all steps, pools and founders, corrects for this reduction.

**Heritability.** Under truncation selection, retaining the top proportion $P_r$ of individuals in replicate $r$ corresponds to a standardised selection intensity $i_r$ (cite Falconer for example). The heritability in the rth replicate, given haplotype frequencies for f founder alleles

$$
h^2_r = \frac{100\,k}{i_r^2} \sum_f \frac{\Delta p_{fr}^2}{p_{fr}}
$$

with $k$ as defined above and $\Delta p_{fr} = \hat{z}_{fr} - \hat{c}_{fr}$.

To obtain an estimate of heritability over all ten replicates we write $e_{fr} = \Delta p_{fr}/i_r$ with mean $\bar{e}_f$ and variance $s^2_f$ 

$$
\hat{h}^2 = 100\,k \sum_f \frac{\bar{e}_f^{\,2} - s_f^2/R}{\hat{p}_f}
$$

where $\hat{p}_f$ is the mean control frequency and $s^2_f/R$ removes the inflation from squaring an estimate. Being unbiased, $\hat{h}^2$ can be negative where the true value is near zero.

**The relationship between the Wald test and $h^2$.** In a case-control genome-wide association study the expected value of the chi-squared statistic is the number of individuals times the variance explained by the causative factor, $E[\chi^2] \approx N \cdot 2p(1-p)\beta^2$. $T$ and $\hat{h}^2$ are likewise both quadratic in the founder frequency shifts, so

$$
E[T] \;\approx\; \frac{n_{\mathrm{eff}}\, \bar{\imath}^{\,2}}{100\,k}\; \hat{h}^2
$$

where $n_{\mathrm{eff}} = \sum_r n_{\mathrm{eff},r}$ and $\bar{\imath}^{\,2}$ is the $n_{\mathrm{eff}}$-weighted mean of $i_r^2$. In Supplementary Figure 1 we plot observed and predicted $\hat{h}^2$ against the Wald statistic and against LWP. At non-significant steps we over-estimate heritability; over moderately significant LWPs of 5 to 15 we over-estimate it by an average of 0.05 percentage points (16-17% relative inflation); and at highly significant steps the estimate is accurate.

The above relationship can also be used to calculate the power of X-QTL experiments. To detect a step contributing 1% of phenotypic variance at an LWP of 15 with 80% power, selecting the longest-lived 5%, one need screen 75,000 flies per treatment. Halving the target heritability to 0.5% doubles the experiment to 150,000 individuals screened.

**Partitioning.** We split the ten replicate contrasts in each treatment into two halves to give eight unbiased estimates of $\hat{h}^2$ at every step — two sexes $\times$ two diets $\times$ two halves. Writing $y_{sdh}$ for the estimate in sex $s$, diet $d$ and half $h$ — a plain symbol, since $\hat{h}^2$ is itself a square — their total sum of squares decomposes into four orthogonal single-degree-of-freedom terms and a four-degree-of-freedom error,

$$
\sum y^2 = 8\bar{y}^2 + SS_{\mathrm{sex}} + SS_{\mathrm{diet}} + SS_{\mathrm{sex:diet}} + SS_{\mathrm{rep}}
$$

where, the design being balanced,

$$
SS_{\mathrm{sex}} = 2(\bar{y}_M - \bar{y}_F)^2, \qquad SS_{\mathrm{diet}} = 2(\bar{y}_{SY20} - \bar{y}_{SY10})^2, \qquad SS_{\mathrm{rep}} = \tfrac{1}{2}\sum_{\mathrm{cells}} (y_{\mathrm{odd}} - y_{\mathrm{even}})^2
$$

$8\bar{y}^2$ is the heritability shared by all four groups, and $SS_{\mathrm{rep}}$ is pure error, the two halves differing only in which replicates they contain. Each term is reported as $MS - MS_{\mathrm{rep}}$, with $MS_{\mathrm{rep}} = SS_{\mathrm{rep}}/4$, and on the heritability scale as $\operatorname{sign}(\cdot)\sqrt{|MS - MS_{\mathrm{rep}}|/8}$, every term being $8\delta^2$ for a deviation $\delta$ on that scale. Since the eight estimates are unbiased so is each contrast, and a term whose true value is near zero falls negative about half the time; a component is taken as present only where its interval excludes zero. The shared term takes its sign from $\bar{y}$, since $8\bar{y}^2$ is a square and a step with strongly negative $\hat{h}^2$ would otherwise read as a large positive shared effect.

The partition is computed at every step, and Figure 1C shows the three components along the genome. To describe patterns genome-wide a coarser unit is needed, since neighbouring steps carry almost the same information — the Wald statistic autocorrelates at $r = 0.999$ over 0.05 cM — and since mapping resolution is genetic rather than physical. The map was therefore cut into *tiles* of 1 cM from the start of each arm, 268 of them covering the 264 cM and 122 Mb of euchromatin. One cM is where neighbouring tiles stop being near-duplicates, adjacent tiles correlating at $r = 0.90$ against $0.97$ at half that size; larger tiles correlate less but merge distinct peaks, 2 cM tiles cutting the separately resolved significant regions from fourteen to eleven. Each tile is summarised by one step, the one within it whose Wald statistic, maximised over the four treatments, is largest, so that all eight estimates refer to the same place; taking the tile mean or its middle step instead changes the partition below by less than the width of its confidence interval. Where nothing has responded the components are all small and their relative sizes arbitrary, so proportions are taken only over the 85 tiles reaching an LWP of 7.5; two of those lie on the X, so the proportions are in effect autosomal. Correlation between tiles falls to $r = 0.38$ only at a separation of 4 cM, so confidence intervals are percentile bootstrap over contiguous groups of four tiles, 2000 resamples.
