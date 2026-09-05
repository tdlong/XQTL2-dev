# What conventional designs could have detected

## The designs

Four ways to map a locus explaining a given fraction of phenotypic variance, in a
base population with `Vp = 1` and `h2 = 50%`, so `Vg = Ve = 0.5`. Autosomes only.

**Outbred panel.** `n` individuals genotyped and phenotyped once, testing a
biallelic SNP at MAF 0.15 on one degree of freedom.

**Inbred lines.** `n` lines, ten measures each, the same one-degree-of-freedom
SNP test on line means.

**Eight-way RILs.** `n` lines, ten measures each, testing eight founder alleles on
seven degrees of freedom.

**X-QTL.** `N` flies screened, the longest-lived proportion `P` selected and
sequenced in pools against a control, testing eight founder haplotypes on seven
degrees of freedom.

## How power is calculated

For the first three, `lambda = n R^2 / (1 - R^2)` with `R^2` the locus
contribution to the variance of one measured unit. For X-QTL the non-centrality
is the relation given in the Methods, `lambda = N P ibar^2 h^2 / (100 k)`. Power
in every case is `P( chi2_df(lambda) > T )`.

For the line-based designs, inbreeding doubles both the locus variance and the
unlinked background, so those cancel; the gain is replication, cutting the
environmental term to `Ve/10`. Hence `Var(line mean) = 2Vg + Ve/10 = 1.05`, and a
locus explaining 2% of individual variance is `R^2 = 0.038` of line variance
against 0.020 for an individual.

Thresholds are LWP 6 for the two SNP panels given the number of tests, LWP 4 for
a RIL linkage scan as established, and LWP 7.5 for X-QTL because that is the
threshold used throughout this paper.

## Power at the sizes these designs are run at

| Design | n | df | Threshold | 2% locus | 1% locus | 0.5% locus |
|:---------------|--------:|---:|:--------|-----:|-----:|------:|
| Outbred panel  |   1,000 |  1 | LWP 6   | 0.35 | 0.04 | 0.004 |
| Inbred lines   |     200 |  1 | LWP 6   | 0.02 | 0.00 | 0.000 |
| Eight-way RILs |     750 |  7 | LWP 4   | 0.70 | 0.16 | 0.022 |
| X-QTL          | 114,000 |  7 | LWP 7.5 | 1.00 | 1.00 | 0.894 |

The 750-line RIL panel is the strongest of the three conventional designs, and at
a 2% locus it has 0.70 power — good, but still missing one such region in three,
on a calculation that assumes a single homogeneous panel.

## The point

X-QTL is the least efficient of these designs per fly, but longevity is phenotyped
in bulk, so it is the only one that can be run at 10^5 individuals — and at the
effect sizes resolved here that is the whole of the difference.
