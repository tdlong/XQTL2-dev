# Supplementary Result X. The power to detect QTL

A comparison of four experimental paradigms for dissecting complex traits in
*Drosophila*, for a locus with heritability $h^2$ in a base population with
phenotypic variance $V_P = 1$ and $V_A = V_E = 0.5$, for autosomes.

**Outbred panel.** $n$ individuals genotyped and phenotyped once, testing a
biallelic SNP at a minor allele frequency of 0.15 on one degree of freedom.

**Inbred lines.** $n$ lines, ten measures per line, the same
one-degree-of-freedom SNP test applied to line means.

**Eight-way RILs.** $n$ lines, ten measures per line, testing eight founder
alleles on seven degrees of freedom.

**X-QTL.** $N$ flies screened, the longest-lived proportion $P$ selected and
sequenced in pools against a contemporaneous control, testing eight founder
haplotypes on seven degrees of freedom.

Under the alternative test statistics are non-central $\chi^2$, so

$$
\mathrm{power} = \Pr\!\left(\chi^2_{\nu}(\lambda) > T_\alpha\right), \qquad
\lambda = \frac{n\,h^2}{1 - h^2}
$$

where $h^2$ is the heritability due to the locus. The line-based designs measure
line means rather than individuals, and take $g^2$ in place of $h^2$. Inbreeding
doubles both the locus variance and the unlinked background, and ten measures per
line cut the environmental term to $V_E/10$, giving

$$
g^{2} = \frac{2\,h^{2}V_P}{2V_A + V_E/10}
$$

so for a locus at $h^2 = 2\%$, $g^2 = 3.8\%$. For X-QTL the non-centrality is the
one given in the Methods. Thresholds are $-\log_{10}P$ of 6 for the two SNP-based
panels reflecting the number of tests carried out, 4 for a RIL linkage scan as
established, and 7.5 for X-QTL. For some common designs the power table below,
and Supplementary Figure 2, illustrate the limited power of traditional
*Drosophila* approaches for dissecting complex traits.

|                |         | h²   |      |       |
|:---------------|--------:|-----:|-----:|------:|
|                | n       | 2%   | 1%   | 0.5%  |
| Outbred panel  |   1,000 | 0.35 | 0.04 | 0.004 |
| Inbred lines   |     200 | 0.02 | 0.00 | 0.000 |
| Eight-way RILs |     750 | 0.70 | 0.16 | 0.022 |
| X-QTL          | 114,000 | 1.00 | 1.00 | 0.894 |

Although X-QTL is the least efficient design per fly, for traits where bulked
phenotyping is an option the overall experiment can be powerful.
