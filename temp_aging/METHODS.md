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

where $n_{\mathrm{eff}} = \sum_r n_{\mathrm{eff},r}$ and $\bar{\imath}^{\,2}$ is the $n_{\mathrm{eff}}$-weighted mean of $i_r^2$. In Supplementary Figure 1 we plot observed and predicted $\hat{h}^2$ against the Wald statistic and against LWP. At non-significant steps we over-estimate heritability; over moderately significant LWPs of 5 to 15 we over-estimate it by an average of 0.05 percentage points (a 16-17% relative inflation); and at highly significant steps the estimate is accurate.

The above relationship also gives the power of an X-QTL experiment. Under the alternative $T$ is non-central $\chi^2_7$ with non-centrality $\lambda = n_{\mathrm{eff}}\bar{\imath}^{\,2}h^2/100k$, so detecting a heritability $h^2$ at threshold $T_\alpha$ with power $1-\beta$ requires $\lambda$ to reach the value $\lambda_{\beta}$ at which $\Pr(\chi^2_7(\lambda) > T_\alpha) = 1-\beta$, that is

$$
n_{\mathrm{eff}} \;=\; \frac{100\,k\,\lambda_{\beta}}{\bar{\imath}^{\,2}\,h^2}, \qquad N \;=\; n_{\mathrm{eff}}/P
$$

flies screened per treatment at selected proportion $P$, since $\bar{\imath}$ is itself fixed by $P$. An LWP of 15 puts $T_\alpha = 85.6$ and $\lambda_{0.2} = 95.6$; selecting the longest-lived 5% gives $\bar{\imath} = 2.06$, so a step contributing 1% of phenotypic variance needs $n_{\mathrm{eff}} = 4{,}500$ and 90,000 flies screened. Requirement scales as $1/h^2$, so halving the target to 0.5% doubles the experiment.

**Partitioning.** We split the ten replicate contrasts in each treatment into two interleaved halves to give eight unbiased estimates of $\hat{h}^2$ at every step — two sexes $\times$ two diets $\times$ two halves. An orthogonal transformation of those eight numbers gives the heritability shared by all four treatments and the parts specific to sex, to diet, and to their interaction, one degree of freedom each, with the difference between the halves of a treatment carrying four. The five terms are orthogonal, so their sums of squares and degrees of freedom add. The interaction is nowhere detectable — it exceeds its 5% critical value at 6% of significant tiles, the rate expected by chance — so it is pooled with the replicate difference, giving an error variance on five degrees of freedom. Writing $\overline{\hat{h}^2}$ for the mean of all eight, $\hat{h}^2_M$, $\hat{h}^2_F$, $\hat{h}^2_{SY10}$, $\hat{h}^2_{SY20}$ for the marginal means, $\hat{h}^2_{(1)}$, $\hat{h}^2_{(2)}$ for the two halves of a treatment and $d_{\mathrm{sex:diet}}$ for the interaction contrast,

$$
d_{\mathrm{shared}} = \overline{\hat{h}^2}, \qquad
d_{\mathrm{sex}} = \tfrac{1}{2}\big(\hat{h}^2_M - \hat{h}^2_F\big), \qquad
d_{\mathrm{diet}} = \tfrac{1}{2}\big(\hat{h}^2_{SY20} - \hat{h}^2_{SY10}\big), \qquad
s^2 = \tfrac{1}{5}\Big[\tfrac{1}{2}\sum_{\mathrm{cells}} \big(\hat{h}^2_{(1)} - \hat{h}^2_{(2)}\big)^2 + 8\,d_{\mathrm{sex:diet}}^2\Big]
$$

Each component is estimated as $\hat{\sigma} = \operatorname{sign}(d)\sqrt{|d^2 - s^2/8|}$ and tested against the error variance as $F_{1,5} = 8d^2/s^2$, with $\operatorname{Var}(d) = s^2/8$ giving the interval $d \pm t_{0.975,5}\sqrt{s^2/8}$. A component below the error is returned negative, which is how a step with nothing there reports.

The partitioning is carried out at every step. To describe genome-wide patterns we also summarize this partitioning and other statistics over 268 non-overlapping one-cM wide tiles spanning the euchromatic genome, with each tile is summarised by the step with the largest Wald statistic maximised over the four treatments. The 1cM tile width is chosen to match the mapping resolution of the synthetic population as two-LWP support intervals on peaks (Table S1) have medians of 0.39 cM above an LWP of 15 and 0.76 cM between 7.5 and 15. Furthermore, 1-cM tiles results in adjacent tiles sharing an $R^2 = 81\%$ of their variance in the Wald statistic. In several cases we summarize statistics over the 85 tiles whose LWP exceeds 7.5, as patterns of variation for non-significant regions are not meaningful.
