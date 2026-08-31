# Methods (incomplete)

Library preparation and sequencing: DNA was prepared on pools of several hundred flies using ?? (Supp table of N of flies).  Illumina DNAseq kit was used for library preparation and was sequenced at UCI Genomics Research and Technology Hub on an Illumina NovaSeq X Plus as 2 × 150 bp reads. Sixty experimental libraries were sequenced: ten replicates of the experiment cages × three groups (unselected control, SY10 selected, SY20 selected) × two sexes.  Reads were aligned to the Drosophila melanogaster dm6 reference genome with bwa mem (Li 2013), coordinate sorted and indexed. Experimental sample alignments were not deduplicated although reads from the eight founder alignments were.

SNP calling: SNPs were called on the five major chromosome arms for the eight founder samples and a catalog of SNPs segregating among the eight founders of the synthetic population created. Founders were called jointly with bcftools (cite): a site entered the catalog if it was biallelic, covered at 10× or more in every founder, had a minor allele frequency of 3% or less in every founder, segregating an alternate allele in at least one founder, and was at least 20 bp from the nearest founder indel. The final catalog comprises 1,207,436 SNPs from 1,887,667 initial candidate positions. For each of the 60 experimental samples reference and alternate read counts were then taken directly from the pileup at each catalog position, and these counts used to infer founder haplotype frequencies. At called SNPs autosomal depth averaged 140× per library and a range of 49× to 400×; fifteen libraries fell below 75× and twenty below 100×. In prior work we have shown 50X coverage is sufficient for high quality haplotype estimation.

Haplotype estimation: The sythetic population is derived from eight inbred founders, so a pool's allele frequency at any SNP is the sum over founders of each founder's frequency in the pool weighted by its genotype. The founder genotypes are known, so over a local window consisting of many SNPs the eight founder frequencies can be estimated via constrained least squares (lsei, R package limSolve) subject to the frequencies summing to one and each being non-negative. Those frequencies were estimated at 5 kb steps, with each step given an initial window of ±75 kb, widened in inverse proportion to the local recombination rate so that windows span comparable genetic rather than physical distance. The fit returns a covariance matrix for the estimates which was averaged over the ten windows either side of each position and used in the Wald test (below). Founders genotypes were hierachically clustered within each window, cut at a height of 2.5, and indistingusihable founders filled by linear interpolation between the means of the flanking resolved windows. The haplotype frequencies were smoothed with a running mean of ±100 kb (41 steps).

**Wald test.** At each window let $\hat{c}$ and $\hat{z}$ be the eight-vectors of estimated founder frequencies in the control and selected pools, and $V_c$, $V_z$ their covariance matrices. The Wald-statistic is

$$
T = \frac{1}{R^2}\,(\hat{c} - \hat{z})^{\top}(V_c + V_z)^{-}(\hat{c} - \hat{z})
$$

on seven degrees of freedom, since the eight haplotype frequencies sum to one. That constraint also makes $V_c + V_z$ singular, so $(\cdot)^{-}$ is a generalised inverse, obtained by eigendecomposition with the null direction discarded and the remaining eigenvalues floored at one per cent of their mean, so that near-degenerate directions cannot inflate $T$. Significance is reported as $-\log_{10} P$.

Each pool's covariance is $V = M + \Sigma$, where $M$ is the multinomial sampling covariance on sampled chromosomes with subscript denoting founder alleles,

$$
M_{ff} = \frac{\bar{p}_f(1-\bar{p}_f)}{n_{\mathrm{eff}}}, \qquad M_{fg} = -\frac{\bar{p}_f \bar{p}_g}{n_{\mathrm{eff}}}
$$

evaluated at the frequency $\bar{p}$ pooled across the selected and control samples, and $\Sigma$ is the smoothed haplotype reconstruction covariance returned by the constrained least squares fit above. The effective number of individuals is obtained estimated from a pool of $n$ flies (where $k$ is one for male X chromosomes and two otherwise) as

$$
n_{\mathrm{eff}} = \frac{kn \sum_f \bar{p}_f(1-\bar{p}_f)}{\sum_f \bar{p}_f(1-\bar{p}_f) + kn \operatorname{tr}(\Sigma)}
$$

The 10 replicates were combined as a weighted using $n_{\mathrm{eff}}$ mean of $\hat{c}$ and $\hat{z}$ with covariances similarly pooled. As the haplotype frequencies entering $T$ are smoothed, their variance is reduced relative unsmoothed estimates. $R^2$ , the squared correlation between raw and smoothed frequencies taken over all windows, pools and founders corrects for this reduction.

**Heritability.** Under truncation selection, retaining the top proportion $P_r$ of individuals in replicate $r$ corresponds to a standardised selection intensity $i_r$ (cite Falconer for example). For a window whose founder haplotype frequencies are $p_f$, that replicate gives

$$
h^2_r = \frac{100\,k}{i_r^2} \sum_f \frac{\Delta p_{fr}^2}{p_f}
$$

with $k$ as defined above and $\Delta p_{fr} = \hat{z}_{fr} - \hat{c}_{fr}$.

To obtain an estimate of heritability over all ten replicates we write $e_{fr} = \Delta p_{fr}/i_r$ with mean $\bar{e}_f$ and variance $s^2_f$ , which over the $R$ replicates gives

$$
\hat{h}^2 = 100\,k \sum_f \frac{\bar{e}_f^{\,2} - s_f^2/R}{\hat{p}_f}
$$

where $\hat{p}_f$ is the mean control frequency and $s^2_f/R$ removes the inflation from squaring an estimate. Being unbiased, $\hat{h}^2$ can be negative where the true value is near zero

OK I am good to here.  And now I have a problem below.  I mean if this is our estimate of h2, this other variance component thing below is out of place?  I think it belong in the variance componet estimation section???

 we also report $h^2_{\mathrm{vc}}$, a variance component fitted per window by maximum likelihood on the scale $\bar{e}_f/\sqrt{\hat{p}_f}$, which is non-negative by construction.

**The relationship between the Wald test and $h^2$.** 

OK, in this section what is the obsession with the non-centrality parameter.  I get that this is what is is.  But for our purposes it is easier to just talk about the relationship between the expected value of the test statistic and underlying causative stuff.  Then the writing is easier to read.  I mean the first equation can be written as the expected value of the X squred statistics.  Then the Wald h2 relationship is just an relationship.  And we never need to belabour the df term.  I mean it is an estimator.  Then you writing is not so forces and defensive.  And it can be short a sweet.  try it.

In a case-control genome-wide association study the non-centrality of the chi-squared test is the number of individuals times the variance explained by the causative factor, that is $N \cdot 2p(1-p)\beta^2$. In the same way, since $T$ and $\hat{h}^2$ are both quadratic in the founder frequency shifts, the non-centrality of the Wald statistic is

$$
T - \mathrm{df} \;\approx\; \frac{n_{\mathrm{eff}}\, \bar{\imath}^{\,2}}{100\,k}\; \hat{h}^2
$$

where $\mathrm{df} = 7$ — a chi-squared statistic averages its degrees of freedom where nothing is happening — $n_{\mathrm{eff}} = \sum_r n_{\mathrm{eff},r}$ and $\bar{\imath}^{\,2}$ is the $n_{\mathrm{eff}}$-weighted mean of $i_r^2$, the same weighting by which the replicates were pooled. The effect reaches the data through the frequency shift truncation selection produces rather than through a measured phenotype, which is why the intensity enters, and a multi-allelic locus contributes $\sum_f \Delta p_f^2/p_f$ in place of $2p(1-p)\beta^2$.

Predicted and observed $\hat{h}^2$ track one another over the range that carries signal (Fig. S1). Between $-\log_{10} P$ of 5 and 15 the observed values run about 0.05 percentage points above the prediction, some 16 to 17% in relative terms, and the two converge above 15.  (This sentence as written as inpentrable ... Just state in Supp Figure 1 we plot predicts h2 against wald statistics as well as out LWP, for non significant regions we over-estimate heritabilty, in the range of LWP of 5 - 15 (significant regions) we over-estimate heritabilty by an average of 0.05 percent (a 16-17% relative inflation), and very accurately estimate h2 for highly signicant regions.



As neff under-estimates kn by a factor of about 1.75 -- this is weird as you expect neff to under-estimate 2N by a factor of 2, so this is confusing.  Like if we are going to do this, why not just state how many individual (n) are required to detect a region contributing 0.5, 1, 2, 4 percent of variation at a LWP of 15 with 80% power?  Like why rearrange A=B/C to AC = B (duh ..) and not really help with intution.

Inverted, it states what a design can resolve: detecting a window of heritability $h^2$ at a threshold $T$ requires

$$
n_{\mathrm{eff}} \approx \frac{100\,k\,(T - \mathrm{df})}{h^2\, \bar{\imath}^{\,2}}
$$

so the detectable effect falls in inverse proportion to the chromosomes sequenced and to the square of the selection intensity. $n_{\mathrm{eff}}$ is below the raw count $kn$ — here by a factor of about 1.75 — because reconstruction error enters the covariance alongside the sampling term, and because the exact weighting is $1/p_f$ rather than $1/[p_f(1-p_f)]$.

**Partitioning.** The ten replicate contrasts per treatment were split into two halves (1 & 2) , giving two independent estimates of heritability at every window in each of the four sex $\times$ diet groups. Writing $y_{sdh}$ for the estimate of heritability for some window in sex $s$, diet $d$ and half $h$, we can express

$$
\sum y^2 = 8\bar{y}^2 + SS_{\mathrm{sex}} + SS_{\mathrm{diet}} + SS_{\mathrm{sex:diet}} + SS_{\mathrm{rep}}
$$

where, the design being balanced,

$$
SS_{\mathrm{sex}} = 2(\bar{y}_M - \bar{y}_F)^2, \qquad SS_{\mathrm{diet}} = 2(\bar{y}_{SY20} - \bar{y}_{SY10})^2, \qquad SS_{\mathrm{rep}} = \tfrac{1}{2}\sum_{\mathrm{cells}} (y_{\mathrm{odd}} - y_{\mathrm{even}})^2
$$

with the treatment terms on one degree of freedom each and $SS_{\mathrm{rep}}$, which is pure error, on four. Each term is reported as $MS - MS_{\mathrm{rep}}$, with $MS_{\mathrm{rep}} = SS_{\mathrm{rep}}/4$, and on the heritability scale as $\operatorname{sign}(MS - MS_{\mathrm{rep}})\sqrt{|MS - MS_{\mathrm{rep}}|/8}$ . Confidence intervals are percentile bootstrap over contiguous groups of tiles.
