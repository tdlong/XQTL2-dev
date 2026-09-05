# What conventional designs could have detected

## The designs

Four ways to map a locus of heritability $h^2$ in a base population with
phenotypic variance $V_P = 1$ and $V_A = V_E = 0.5$. Autosomes only.

**Outbred panel.** $n$ individuals genotyped and phenotyped once, testing a
biallelic SNP at a minor allele frequency of 0.15 on one degree of freedom.

**Inbred lines.** $n$ lines, ten measures each, the same one-degree-of-freedom
SNP test applied to line means.

**Eight-way RILs.** $n$ lines, ten measures each, testing eight founder alleles on
seven degrees of freedom.

**X-QTL.** $N$ flies screened, the longest-lived proportion $P$ selected and
sequenced in pools against a contemporaneous control, testing eight founder
haplotypes on seven degrees of freedom.

## How power is calculated

Under the alternative every test statistic is non-central $\chi^2$, so

$$
\mathrm{power} \;=\; \Pr\!\left(\chi^2_{\nu}(\lambda) > T_\alpha\right)
$$

$h^2$ is the heritability due to the locus, the variance explained by the
causative factor. For the outbred panel and for X-QTL it enters directly. The two
line-based designs measure line means instead: inbreeding doubles both the locus
variance and the unlinked background, so those cancel, and ten measures per line
cut the environmental term to $V_E/10$. Write

$$
{}^{\prime}h^{2} \;=\; \frac{2\,h^{2}V_P}{2V_A + V_E/10}
$$

for the locus on that scale, so a locus at $h^2 = 2\%$ is
${}^{\prime}h^{2} = 3.8\%$. The non-centrality is then

$$
\lambda = \frac{n\,h^{2}}{1-h^{2}}, \qquad
\lambda = \frac{n\,{}^{\prime}h^{2}}{1-{}^{\prime}h^{2}}, \qquad
\lambda = \frac{N P \bar{\imath}^{\,2} h^{2}}{100\,k}
$$

for the outbred panel, for the two line-based designs, and for X-QTL as given in
the Methods.

Thresholds are $-\log_{10}P$ of 6 for the two SNP panels given the number of tests,
4 for a RIL linkage scan as established, and 7.5 for X-QTL, that being the
threshold used throughout this paper.

## Power at the sizes these designs are run at

| Design | $n$ | $\nu$ | $h^2 = 2\%$ | $h^2 = 1\%$ | $h^2 = 0.5\%$ |
|:---------------|--------:|---:|-----:|-----:|------:|
| Outbred panel  |   1,000 |  1 | 0.35 | 0.04 | 0.004 |
| Inbred lines   |     200 |  1 | 0.02 | 0.00 | 0.000 |
| Eight-way RILs |     750 |  7 | 0.70 | 0.16 | 0.022 |
| X-QTL          | 114,000 |  7 | 1.00 | 1.00 | 0.894 |

The 750-line RIL panel is the strongest of the three conventional designs, and at
$h^2 = 2\%$ it has 0.70 power — good, but still missing one such region in three,
on a calculation that assumes a single homogeneous panel.

## The point

X-QTL is the least efficient of these designs per fly, but longevity is phenotyped
in bulk, so it is the only one that can be run at $10^5$ individuals — and at the
effect sizes resolved here that is the whole of the difference.
