# Methods (incomplete)

Library preparation and sequencing: DNA was prepared on pools of several hundred flies using ?? (Supp table of N of flies).  Illumina DNAseq kit was used for library preparation and was sequenced at UCI Genomics Research and Technology Hub on an Illumina NovaSeq X Plus as 2 × 150 bp reads. Sixty experimental libraries were sequenced: ten replicates of the experiment cages × three groups (unselected control, SY10 selected, SY20 selected) × two sexes.  Reads were aligned to the Drosophila melanogaster dm6 reference genome with bwa mem (Li 2013), coordinate sorted and indexed. Experimental sample alignments were not deduplicated although reads from the eight founder alignments were.

SNP calling: SNPs were called on the five major chromosome arms for the eight founder samples and a catalog of SNPs segregating among the eight founders of the synthetic population created. Founders were called jointly with bcftools (cite): a site entered the catalog if it was biallelic, covered at 10× or more in every founder, had a minor allele frequency of 3% or less in every founder, segregating an alternate allele in at least one founder, and was at least 20 bp from the nearest founder indel. The final catalog comprises 1,207,436 SNPs from 1,887,667 initial candidate positions. For each of the 60 experimental samples reference and alternate read counts were then taken directly from the pileup at each catalog position, and these counts used to infer founder haplotype frequencies. At called SNPs autosomal depth averaged 140× per library and a range of 49× to 400×; fifteen libraries fell below 75× and twenty below 100×. In prior work we have shown 50X coverage is sufficient for high quality haplotype estimation.

Haplotype estimation: The sythetic population is derived from eight inbred founders, so a pool's allele frequency at any SNP is the sum over founders of each founder's frequency in the pool weighted by its genotype. The founder genotypes are known, so over a local window consisting of many SNPs the eight founder frequencies can be estimated via constrained least squares (lsei, R package limSolve) subject to the frequencies summing to one and each being non-negative. Those frequencies were estimated at 5 kb steps, with each step given an initial window of ±75 kb, widened in inverse proportion to the local recombination rate so that windows span comparable genetic rather than physical distance. The fit returns a covariance matrix for the estimates which was averaged over the ten windows either side of each position and used in the Wald test (below). Founders genotypes were hierachically clustered within each window, cut at a height of 2.5, and indistingusihable founders filled by linear interpolation between the means of the flanking resolved windows. The haplotype frequencies were smoothed with a running mean of ±100 kb (41 steps).

**Wald test.** At each window let $\hat{c}$ and $\hat{z}$ be the eight-vectors of estimated founder frequencies in the control and selected pools, and $V_c$, $V_z$ their covariance matrices. The Wald-statistic is

$$
T = (\hat{c} - \hat{z})^{\top}\,(V_c + V_z)^{-1}\,(\hat{c} - \hat{z})
$$

distributed under the null as $\chi^2$ on seven degrees of freedom (and $V_c + V_z$ is correspondingly singular. The null eigenvector was dropped and the remaining eigenvalues floored at one per cent of their mean before inversion -- I am confused here, we seem to refer to eigenvectors before inversion of var/cov matrices, but where do these come from without some sort of SVD?? - so this is not clear -- I think the 7 degrees of freedom is obvious).

Each pool's covariance is $V = M + \Sigma$, where $M$ is the multinomial sampling convariance on sampled chromosomes,

$$
M_{ff} = \frac{\bar{p}_f(1-\bar{p}_f)}{n_{\mathrm{eff}}}, \qquad M_{fg} = -\frac{\bar{p}_f \bar{p}_g}{n_{\mathrm{eff}}}
$$

evaluated at the frequency $\bar{p}$ pooled across the selected and control samples, and $\Sigma$ is the haplotype reconstruction covariance returned by the constrained least squares fit above, averaged over the ten windows either side. The effective number of sampled chromosomes is

$$
n_{\mathrm{eff}} = \frac{nS}{S + nC}, \qquad S = \sum_f \bar{p}_f(1-\bar{p}_f), \qquad C = \operatorname{tr}(\Sigma)
$$

with $n$ the number of chromosomes sampled: two per fly on an autosome or the X in females, one per fly on the hemizygous X in males.   (I do not love this for a few reasons, why not just use the kn notation and define n as the number of individuals in the pool and k as one for male X chromsomes, and 2 otherwise.  I further am not sure I love the choice of symbols S and C?  As you might imagine S and C are related to selected and control individuals... whereas they clearly represent other things conceptuallly, so perhaps we could discuss better nomenclature).  Replicate cages were combined as a weighted mean of $\hat{c}$ and $\hat{z}$ with weights $n_{\mathrm{eff}}$, their covariances pooled to match. $T$ was divided by the mean squared correlation between raw and smoothed founder frequencies to account for the smoothing, and significance is reported as $-\log_{10} P$.  (not loving this either ... as we sort of define T right at the start of the paragraph ... then in the last sentence we pull a fast one and add a normalizing constant, so I think the initial definition of T needs like a rho term, and then we define rho as the end as a term to account for the fact we using smoothed founder haplotype frequencies.)

**Heritability.** Under truncation selection for any given replicate retaining the top proportion $P$ of individuals corresponds to a standardised selection intensity $i$ (cite Falconer for example), and for some window a founder haplotype at frequency $p_f$ with additive effect $a_f$ has a difference in frequency between selected and control pools of $\Delta p_f = p_f (a_f - \bar{a}) i$, where $\bar{a}$ is the mean effect, so that $a_f - \bar{a} = \Delta p_f / (p_f i)$ (this is a strange way to do this, why not just where a bar is and give the equation, like why express it in this form?????).

An individual's genotypic value at the locus is the sum over the allele copies it carries, so the additive variance contributed by the window is the per-copy variance $\sum_f p_f (a_f - \bar{a})^2$ multiplied by the number of copies. On an autosome a fly carries two, and substituting the expression above,

$$
V_A = \frac{100k}{Vt i^2} \sum_f \frac{\Delta p_f^2}{p_f}
$$

(then you lost me here.  I mean once we define delta p f, we can just say and the additive variance is defined as... ).  Then the next sentence is awkward, when we could just say given Va the heritability for the window is Va/Vt expressed as a percentage.



The trait being standardised, $H^2 = V_A$, reported as a percentage.

Then you go off the rails.  Maybe we don;t even need to define Va, lets just define heritability in one go, and instead of that "2" use a k (and we can note as defined above) and add the Vt and we are done?  Like it seems we can just have one simpel equation.  I even gave it a go in the equation although you have to fix up the formatting).



A male carries one X rather than two, so at an X-linked locus in males the same set of frequency shifts corresponds to half the additive variance and the leading constant is 100 rather than 200. Females carry two X chromosomes and take the autosomal form. Writing $k$ for the number of copies, the constant is $100k$ throughout.

Then we can say how we average over the replicates.  I mean we went to all this trouble to define Va ( or h2) ... then you don't use it....  so that is sort of awkward in terms of presentations.???  And I think we need to use a lower case H2.

Replicates differ in the proportion selected, so the quantity they measure in common is the response per unit selection intensity, $e_{fr} = (\hat{z}_{fr} - \hat{c}_{fr})/i_r$, of which the $n$ replicates are independent estimates. Their mean $\bar{e}_f$ estimates $p_f(a_f - \bar{a})$, and its square is biased upward by the sampling variance of the mean, $s_f^2/n$, giving

$$
\hat{H}^2 = 100\,k \sum_f \frac{\bar{e}_f^{\,2} - s_f^2/n}{\hat{p}_f}
$$

with $\hat{p}_f$ the mean control frequency and $s_f^2$ the variance of $e_{fr}$ across replicates. The bias correction is measured from the replicates, so the estimator requires no model of the sampling or reconstruction covariances. Being unbiased it takes negative values where the true value is near zero; we therefore also report $H^2_{\mathrm{vc}}$, in which a single variance component is fitted per window by maximum likelihood on the scale $\bar{e}_f/\sqrt{\hat{p}_f}$ with known per-founder noise, so that non-negativity is a boundary solution of the likelihood rather than a truncation.

**Relation between the test and the effect size.** The Wald statistic and $H^2$ are the same quadratic form in the founder frequency shifts and differ only in their weighting. Where the covariance is dominated by multinomial sampling, $V \approx p_f(1-p_f)/n$, so $T \approx n \sum_f \Delta p_f^2 / [p_f(1-p_f)]$, against $H^2 = (100k/i^2)\sum_f \Delta p_f^2/p_f$. The two are therefore proportional,

$$
H^2 \approx \frac{100\,k}{n_{\mathrm{eff}}\, i^2}\; T
$$

with the constant fixed within a scan by the number of chromosomes sampled and the intensity of selection, and nothing else. This expression is analgous to the familiar non-centrality parameter from case control GWAS designs where $\chi^2 \approx N \cdot 2p(1-p)\beta^2$ . (We may have to further clarify ... as there are not measured phenotypes in the case control, so it is unclear).  The difference is that the effect reaches the data through truncation selection rather than through measured phenotypes, so the selection intensity enters as $i^2$, and that a multi-allelic locus contributes $\sum_f \Delta p_f^2/p_f$ rather than a single $2p(1-p)\beta^2$.

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
