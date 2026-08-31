# Methods (incomplete)

Library preparation and sequencing: DNA was prepared on pools of several hundred flies using ?? (Supp table of N of flies).  Illumina DNAseq kit was used for library preparation and was sequenced at UCI Genomics Research and Technology Hub on an Illumina NovaSeq X Plus as 2 × 150 bp reads. Sixty experimental libraries were sequenced: ten replicates of the experiment cages × three groups (unselected control, SY10 selected, SY20 selected) × two sexes.  Reads were aligned to the Drosophila melanogaster dm6 reference genome with bwa mem (Li 2013), coordinate sorted and indexed. Experimental sample alignments were not deduplicated although reads from the eight founder alignments were.

SNP calling: SNPs were called on the five major chromosome arms for the eight founder samples and a catalog of SNPs segregating among the eight founders of the synthetic population created. Founders were called jointly with bcftools (cite): a site entered the catalog if it was biallelic, covered at 10× or more in every founder, had a minor allele frequency of 3% or less in every founder, segregating an alternate allele in at least one founder, and was at least 20 bp from the nearest founder indel. The final catalog comprises 1,207,436 SNPs from 1,887,667 initial candidate positions. For each of the 60 experimental samples reference and alternate read counts were then taken directly from the pileup at each catalog position, and these counts used to infer founder haplotype frequencies. At called SNPs autosomal depth averaged 140× per library and a range of 49× to 400×; fifteen libraries fell below 75× and twenty below 100×. In prior work we have shown 50X coverage is sufficient for high quality haplotype estimation.

Haplotype estimation: The sythetic population is derived from eight inbred founders, so a pool's allele frequency at any SNP is the sum over founders of each founder's frequency in the pool weighted by its genotype. The founder genotypes are known, so over a local window consisting of many SNPs the eight founder frequencies can be estimated via constrained least squares (lsei, R package limSolve) subject to the frequencies summing to one and each being non-negative. Those frequencies were estimated at 5 kb steps, with each step given an initial window of ±75 kb, widened in inverse proportion to the local recombination rate so that windows span comparable genetic rather than physical distance. The fit returns a covariance matrix for the estimates which was averaged over the ten windows either side of each position and used in the Wald test (below). Founders genotypes were hierachically clustered within each window, cut at a height of 2.5, and indistingusihable founders filled by linear interpolation between the means of the flanking resolved windows. The haplotype frequencies were smoothed with a running mean of ±100 kb (41 steps).

**Wald test.** At each window let $\hat{c}$ and $\hat{z}$ be the eight-vectors of estimated founder frequencies in the control and selected pools, and $V_c$, $V_z$ their covariance matrices. The Wald statistic is

$$
T = \frac{1}{\rho}\,(\hat{c} - \hat{z})^{\top}(V_c + V_z)^{-}(\hat{c} - \hat{z})
$$

where $(\cdot)^{-}$ is a generalised inverse and $\rho$ corrects for smoothing, both defined below. Because the eight haplotype frequencies sum to one, $V_c + V_z$ is singular and $T$ carries seven degrees of freedom. The generalized inverse is taken by eigendecomposition, discarding the null direction and flooring the remaining eigenvalues at one per cent of their mean so that near-degenerate directions cannot inflate $T$. Under the null $T$ is $\chi^2_7$, and significance is reported as $-\log_{10} P$.

Each pool's covariance is $V = M + \Sigma$. $M$ is the multinomial covariance of the chromosomes sampled,

$$
M_{ff} = \frac{\bar{p}_f(1-\bar{p}_f)}{n_{\mathrm{eff}}}, \qquad M_{fg} = -\frac{\bar{p}_f\bar{p}_g}{n_{\mathrm{eff}}}
$$

evaluated at the frequency $\bar{p}$ pooled across the selected and control samples, and $\Sigma$ is the haplotype reconstruction covariance returned by the constrained least squares fit above, averaged over the ten windows either side. A pool of $n$ flies carrying $k$ copies of the locus each supplies $kn$ chromosomes (where $k$ is one for male X chromsomes and two otherwise), which reconstruction error reduces to an effective count

$$
n_{\mathrm{eff}} = \frac{kn \sum_f \bar{p}_f(1-\bar{p}_f)}{\sum_f \bar{p}_f(1-\bar{p}_f) + kn\operatorname{tr}(\Sigma)}
$$

Replicate cages were combined as a weighted mean of $\hat{c}$ and $\hat{z}$ with weights $n_{\mathrm{eff}}$, their covariances pooled to match. Our smoothing of the founder frequencies removes variance from them, and $\rho$, the mean squared correlation between the raw and smoothed frequencies, restores it.

**Heritability.** For a replicate retaining the longest-lived proportion $P$ of individuals, the standardised selection intensity is $i$ (Falconer and Mackay). At a window, a founder haplotype at frequency $p_f$ whose additive effect on the trait is $a_f$, expressed in phenotypic standard deviations, differs in frequency between the selected and control pools by

$$
\Delta p_f = p_f (a_f - \bar{a})\, i
$$

with $\bar{a}$ the mean effect over founders. An individual's genotypic value at the locus is the sum over the $k$ copies it carries, so the window contributes additive variance $k\sum_f p_f (a_f - \bar{a})^2$, and because $a_f$ is in phenotypic standard deviations this is already a fraction of the phenotypic variance. Substituting $a_f - \bar{a} = \Delta p_f/(p_f i)$, the heritability of the window as a percentage is

$$
h^2 = \frac{100\,k}{i^2} \sum_f \frac{\Delta p_f^2}{p_f}
$$

$\Delta p_f / i = p_f(a_f - \bar{a})$ is a property of the window and not of the replicate, so every replicate estimates the same quantity however stringently it was selected. Writing $e_{fr} = (\hat{z}_{fr} - \hat{c}_{fr})/i_r$ for replicate $r$, and $\bar{e}_f$ and $s^2_f$ for the mean and variance of $e_{fr}$ over the $R$ replicates, the expression above becomes $h^2 = 100k\sum_f \bar{e}_f^{\,2}/\hat{p}_f$. Squaring an estimate inflates it by its own sampling variance, $s^2_f/R$, and subtracting that gives

$$
\hat{h}^2 = 100\,k \sum_f \frac{\bar{e}_f^{\,2} - s_f^2/R}{\hat{p}_f}
$$

with $\hat{p}_f$ the mean control frequency. The correction is measured from the replicates, so the estimator needs no model of the sampling or reconstruction covariances. Being unbiased it takes negative values where the true value is near zero; we therefore also report $h^2_{\mathrm{vc}}$, in which a single variance component is fitted per window by maximum likelihood on the scale $\bar{e}_f/\sqrt{\hat{p}_f}$ with known per-founder noise, so that non-negativity is a boundary solution of the likelihood rather than a truncation.

**The relationship between the Wald test and $h^2$.** $T$ and $h^2$ are both quadratic in the founder frequency shifts and differ only in their weighting: where the covariance is dominated by multinomial sampling $T \approx n_{\mathrm{eff}}\sum_f \Delta p_f^2/[p_f(1-p_f)]$, against $h^2 = (100k/i^2)\sum_f \Delta p_f^2/p_f$. In a case-control genome-wide association study the same argument gives the non-centrality parameter $\chi^2 \approx N\cdot 2p(1-p)\beta^2$, the number of individuals times the variance the locus explains. Arranged the same way, this design gives

$$
T \approx \frac{n_{\mathrm{eff}}\, \bar{\imath}^{\,2}}{100\,k}\; h^2
$$

with $n_{\mathrm{eff}}$ the effective chromosome count pooled over replicates and $\bar{\imath}$ their mean selection intensity. The effect reaches the data through the frequency shift truncation selection produces rather than through a measured phenotype, which is why the intensity enters, and a multi-allelic locus contributes $\sum_f \Delta p_f^2/p_f$ in place of $2p(1-p)\beta^2$.

The relation holds in the data: over windows with $-\log_{10} P > 5$ the ratio $h^2/T$ is 0.0077 in autosomal females, 0.0076 in autosomal males and 0.0072 on the male X, one constant across arms and sexes including the male X, where $k$ and $n_{\mathrm{eff}}$ both halve and cancel.

Inverted, it states what a design can resolve: detecting a window of heritability $h^2$ at a threshold $T$ requires

$$
n_{\mathrm{eff}} \approx \frac{100\,k\,T}{h^2 \bar{\imath}^{\,2}}
$$

so the detectable effect falls in inverse proportion to the chromosomes sequenced and to the square of the selection intensity. $n_{\mathrm{eff}}$ is below the raw count $kn$ — here by a factor of about 1.75 — because reconstruction error enters the covariance alongside the sampling term, and because the exact weighting is $1/p_f$ rather than $1/[p_f(1-p_f)]$. For the present design $-\log_{10} P = 7.5$ corresponds to 0.38% of phenotypic variance.

**Partitioning.** The ten replicate contrasts per treatment were split into odd and even sets of five, giving two independent estimates of heritability at every window in each of the four sex $\times$ diet groups. Writing $y_{sdh}$ for the estimate in sex $s$, diet $d$ and half $h$, the eight values were decomposed as a balanced $2 \times 2$ with two replicates per cell, uncorrected for the mean,

$$
\sum y^2 = 8\bar{y}^2 + SS_{\mathrm{sex}} + SS_{\mathrm{diet}} + SS_{\mathrm{sex:diet}} + SS_{\mathrm{rep}}
$$

where, the design being balanced,

$$
SS_{\mathrm{sex}} = 2(\bar{y}_M - \bar{y}_F)^2, \qquad SS_{\mathrm{diet}} = 2(\bar{y}_{SY20} - \bar{y}_{SY10})^2, \qquad SS_{\mathrm{rep}} = \tfrac{1}{2}\sum_{\mathrm{cells}} (y_{\mathrm{odd}} - y_{\mathrm{even}})^2
$$

with the treatment terms on one degree of freedom each and $SS_{\mathrm{rep}}$, which is pure error, on four. Each term is reported as $MS - MS_{\mathrm{rep}}$, with $MS_{\mathrm{rep}} = SS_{\mathrm{rep}}/4$, and on the heritability scale as $\operatorname{sign}(MS - MS_{\mathrm{rep}})\sqrt{|MS - MS_{\mathrm{rep}}|/8}$. Terms were not truncated at zero. Confidence intervals are percentile bootstrap over contiguous groups of tiles.
