# Methods (incomplete)

Library preparation and sequencing: DNA was prepared on bulked flies using ?? (table of N of flies).  Library preparation used Illumina DNAseq kit and was sequenced at UCI Genomics Research and Technology Hub on an Illumina NovaSeq X Plus as 2 × 150 bp reads. Sixty experimental libraries were sequenced: ten replicates of the experiment cages × three groups (unselected control, SY10 selected, SY20 selected) × two sexes.  Reads were aligned to the Drosophila melanogaster dm6 reference genome with bwa mem (Li 2013), coordinate sorted and indexed. Experimental sample alignments were not deduplicated although reads from the eight founder alignments were.

SNP calling: SNPs were called on the five major chromosome arms for the eight founder samples and a catalog of SNPs segregating among the eight founders of the synthetic population created. Founders were called jointly with bcftools, and a site entered the catalog if it was biallelic; covered at 10× or more in every founder; had a minor allele frequency of 3% or less in every founder; segregating an alternate allele in at least one founder; and at least 20 bp from the nearest founder indel. Founder B5 was exempted from these criteria on chromosome 2L (as the bam is created from inversion free alleles). The final catalog comprises 1,207,436 SNPs from 1,887,667 initial candidate positions. For each of the 60 experimental sample reference and alternate read counts for each pooled library were then taken directly from the pileup at each catalog position, and these counts used to infer founder haplotype frequencies. At called SNPs autosomal depth averaged 140× per library and a range of 49× to 400×; fifteen libraries fell below 75× and twenty below 100×. In prior work we have shown 50X coverage is sufficient for high quality haplotype estimation.

Haplotype estimation: The sythetic population is derived from eight inbred founders, so a pool's allele frequency at any SNP is a mixture: the sum over founders of each founder's frequency in the pool, weighted by its genotype at that SNP. The founder genotypes are known, so with many SNPs in a local window the eight founder frequencies can be recovered by least squares. Founder frequencies were estimated at 5 kb steps. Each step was given an initial window of ±75 kb, widened in inverse proportion to the local recombination rate so that windows span comparable genetic rather than physical distance; recombination rate was taken from an eighth degree polynomial fitted to each chromosome arm. Windows holding fewer than 50 catalog SNPs were dropped. At each window, and for each library independently, the eight frequencies were estimated by constrained least squares (lsei, R package limSolve), minimising the squared deviation between observed and predicted allele frequency subject to the frequencies summing to one and each being non-negative. The fit also returns a covariance matrix for the estimates; this covariance matrix was averaged over the ten windows either side of each position and carried into the tests. Founders cannot always be told apart. Where two or more carry the same haplotype across a window their individual frequencies are unidentifiable and only their sum is determined. Founders were therefore clustered within each window by hierarchical clustering of their genotypes, cut at a height of 2.5, and any founder sharing a cluster with another was masked. Masked windows, together with the few in which the fit failed to satisfy its own constraints, were filled by linear interpolation between the means of the flanking resolved windows. The filled series was then smoothed with a running mean of ±100 kb (41 steps).

**Wald test.** At each window let $\hat{c}$ and $\hat{z}$ be the eight-vectors of founder frequencies in the control and selected pools, and $V_c$, $V_z$ their covariance matrices. The statistic is

$$
T = (\hat{c} - \hat{z})^{\top}\,(V_c + V_z)^{-1}\,(\hat{c} - \hat{z})
$$

distributed under the null as $\chi^2$ on seven degrees of freedom, one fewer than the number of founders because the frequencies sum to one and $V_c + V_z$ is correspondingly singular. The null eigenvector was dropped and the remaining eigenvalues floored at one per cent of their mean before inversion.

Each pool's covariance is $V = M + \Sigma$, where $M$ is the multinomial sampling covariance,

$$
M_{ff} = \frac{\bar{p}_f(1-\bar{p}_f)}{n_{\mathrm{eff}}}, \qquad M_{fg} = -\frac{\bar{p}_f \bar{p}_g}{n_{\mathrm{eff}}}
$$

evaluated at the frequency $\bar{p}$ pooled across the selected and control samples, and $\Sigma$ is the haplotype reconstruction covariance returned by the constrained least squares fit, averaged over the ten windows either side. The effective number of chromosomes is

$$
n_{\mathrm{eff}} = \frac{nS}{S + nC}, \qquad S = \sum_f \bar{p}_f(1-\bar{p}_f), \qquad C = \operatorname{tr}(\Sigma)
$$

with $n$ the number of chromosomes sampled: two per fly on an autosome or the X in females, one per fly on the hemizygous X in males. Replicate cages were combined as a weighted mean of $\hat{c}$ and $\hat{z}$ with weights $n_{\mathrm{eff}}$, their covariances pooled to match. $T$ was divided by the mean squared correlation between raw and smoothed founder frequencies to account for the smoothing, and significance is reported as $-\log_{10} P$.

**Heritability.** Under truncation selection retaining the top proportion $P$ with standardised selection intensity $i$, a founder haplotype at frequency $p_f$ with additive effect $a_f$ changes in frequency by $\Delta p_f = p_f (a_f - \bar{a}) i$, where $\bar{a}$ is the mean effect, so that $a_f - \bar{a} = \Delta p_f / (p_f i)$.

An individual's genotypic value at the locus is the sum over the allele copies it carries, so the additive variance contributed by the window is the per-copy variance $\sum_f p_f (a_f - \bar{a})^2$ multiplied by the number of copies. On an autosome a fly carries two, and substituting the expression above,

$$
V_A = \frac{2}{i^2} \sum_f \frac{\Delta p_f^2}{p_f}
$$

The trait being standardised, $H^2 = V_A$, reported as a percentage.

A male carries one X rather than two, so at an X-linked locus in males the same set of frequency shifts corresponds to half the additive variance and the leading constant is 100 rather than 200. Females carry two X chromosomes and take the autosomal form. Writing $k$ for the number of copies, the constant is $100k$ throughout.

Replicates differ in the proportion selected, so the quantity they measure in common is the response per unit selection intensity, $e_{fr} = (\hat{z}_{fr} - \hat{c}_{fr})/i_r$, of which the $n$ replicates are independent estimates. Their mean $\bar{e}_f$ estimates $p_f(a_f - \bar{a})$, and its square is biased upward by the sampling variance of the mean, $s_f^2/n$, giving

$$
\hat{H}^2 = 100\,k \sum_f \frac{\bar{e}_f^{\,2} - s_f^2/n}{\hat{p}_f}
$$

with $\hat{p}_f$ the mean control frequency and $s_f^2$ the variance of $e_{fr}$ across replicates. The bias correction is measured from the replicates, so the estimator requires no model of the sampling or reconstruction covariances. Being unbiased it takes negative values where the true value is near zero; we therefore also report $H^2_{\mathrm{vc}}$, in which a single variance component is fitted per window by maximum likelihood on the scale $\bar{e}_f/\sqrt{\hat{p}_f}$ with known per-founder noise, so that non-negativity is a boundary solution of the likelihood rather than a truncation.

**Relation between the test and the effect size.** The Wald statistic and $H^2$ are the same quadratic form in the founder frequency shifts and differ only in their weighting. Where the covariance is dominated by multinomial sampling, $V \approx p_f(1-p_f)/n$, so $T \approx n \sum_f \Delta p_f^2 / [p_f(1-p_f)]$, against $H^2 = (100k/i^2)\sum_f \Delta p_f^2/p_f$. The two are therefore proportional,

$$
H^2 \approx \frac{100\,k}{n_{\mathrm{eff}}\, i^2}\; T
$$

with the constant fixed within a scan by the number of chromosomes sampled and the intensity of selection, and nothing else. This is the non-centrality parameter for the design: as in a genome-wide association study, where $\chi^2 \approx N \cdot 2p(1-p)\beta^2$ ties the test statistic to sample size and variance explained, the statistic here is the product of the number of chromosomes sequenced and the heritability of the window. The difference is that the effect reaches the data through truncation selection rather than through measured phenotypes, so the selection intensity enters as $i^2$, and that a multi-allelic locus contributes $\sum_f \Delta p_f^2/p_f$ rather than a single $2p(1-p)\beta^2$.

This holds in the data. Over windows with $-\log_{10} P > 5$, the ratio $H^2/T$ is 0.0077 in autosomal females, 0.0076 in autosomal males and 0.0072 on the male X — one constant across arms and sexes, the male X included, where both $k$ and $n$ are halved and cancel. Taken the other way it predicts the heritability at a given significance threshold: $T = 48.3$ at $-\log_{10} P = 7.5$ on seven degrees of freedom gives 0.37%, against 0.38% observed, and $T = 85.6$ at $-\log_{10} P = 15$ gives 0.65% against 0.68%.

Inverted, it states what a design can resolve. Detecting a window of heritability $H^2$ at a threshold $T$ requires

$$
n_{\mathrm{eff}} \approx \frac{100\,k\,T}{H^2 i^2}
$$

so the detectable effect falls in inverse proportion to the number of chromosomes sequenced and in inverse proportion to the square of the selection intensity. $n_{\mathrm{eff}}$ is smaller than the raw chromosome count — here by a factor of about 1.75 — because the haplotype reconstruction error enters the covariance alongside the sampling term, and because the exact weighting is $1/p_f$ rather than $1/[p_f(1-p_f)]$. For the present design the threshold of $-\log_{10} P = 7.5$ corresponds to 0.38% of phenotypic variance.

**Partitioning.** The ten replicate contrasts per treatment were split into odd and even sets of five, giving two independent estimates of heritability at every window in each of the four sex $\times$ diet groups. Writing $y_{sdh}$ for the estimate in sex $s$, diet $d$ and half $h$, the eight values were decomposed as a balanced $2 \times 2$ with two replicates per cell, uncorrected for the mean,

$$
\sum y^2 = 8\bar{y}^2 + SS_{\mathrm{sex}} + SS_{\mathrm{diet}} + SS_{\mathrm{sex:diet}} + SS_{\mathrm{rep}}
$$

where, the design being balanced,

$$
SS_{\mathrm{sex}} = 2(\bar{y}_M - \bar{y}_F)^2, \qquad SS_{\mathrm{diet}} = 2(\bar{y}_{SY20} - \bar{y}_{SY10})^2, \qquad SS_{\mathrm{rep}} = \tfrac{1}{2}\sum_{\mathrm{cells}} (y_{\mathrm{odd}} - y_{\mathrm{even}})^2
$$

with the treatment terms on one degree of freedom each and $SS_{\mathrm{rep}}$, which is pure error, on four. Each term is reported as $MS - MS_{\mathrm{rep}}$, with $MS_{\mathrm{rep}} = SS_{\mathrm{rep}}/4$, and on the heritability scale as $\operatorname{sign}(MS - MS_{\mathrm{rep}})\sqrt{|MS - MS_{\mathrm{rep}}|/8}$. Terms were not truncated at zero. Confidence intervals are percentile bootstrap over contiguous groups of tiles.
