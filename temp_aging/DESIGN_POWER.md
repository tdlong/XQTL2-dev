# What conventional designs could have detected

Autosomes only. Computed by `design_power.R`; the figure is `make_FigS2.R`.

## The designs

Each design detects a locus explaining a fraction of phenotypic variance in a base
population with `Vp = 1` and `h2 = 50%`, so `Vg = Ve = 0.5`. Power is

    power = P( chi2_df(lambda) > T )      lambda = n R^2 / (1 - R^2)

with `n` the number of units measured and `R^2` the locus contribution to the
variance of one such unit. Inbreeding doubles both the locus variance and the
unlinked background, so those cancel; the gain from lines is replication, which
cuts the environmental term to `Ve/10`. Hence `Var(line mean) = 2Vg + Ve/10 = 1.05`
and a locus explaining 2% of individual variance is `R^2 = 0.038` of line variance.

| Design | Unit | df | Threshold | `R^2` at a 2% locus |
|:-------------------|:-------------------------|---:|:-----------|--------:|
| Outbred panel      | individual, one measure  |  1 | LWP 6      |  0.0200 |
| Inbred lines       | line, ten measures       |  1 | LWP 6      |  0.0381 |
| Eight-way RILs     | line, ten measures       |  7 | LWP 4      |  0.0381 |
| X-QTL              | fly screened             |  7 | LWP 7.5    |       — |

The two SNP designs test a biallelic locus at MAF 0.15 on one degree of freedom;
the eight-way designs test eight founder alleles on seven. Thresholds are LWP 6
for the SNP panels given the number of tests, LWP 4 for a RIL linkage scan as
established, and LWP 7.5 for X-QTL because that is the threshold used throughout
this paper. For X-QTL the non-centrality is the relation given in the Methods,
`lambda = N P ibar^2 h^2 / (100 k)`, with `N` flies screened at selected
proportion `P`.

## Power at the sizes these designs are run at

| Design | n | Power, 2% locus | Power, 1% locus | Power, 0.5% locus |
|:-------------------|--------:|-----:|-----:|------:|
| Outbred panel      |   1,000 | 0.35 | 0.04 | 0.004 |
| Inbred lines       |     200 | 0.02 | 0.00 | 0.000 |
| Eight-way RILs     |     750 | 0.70 | 0.16 | 0.022 |
| X-QTL              | 114,000 | 1.00 | 1.00 | 0.894 |

## Supplementary Figure 2

Power against experiment size, at a locus explaining 2% of phenotypic variance
(left) and 0.5%, the size resolved here (right). The x axis is individuals
phenotyped, or lines for the inbred-line and RIL panels. Solid segments span the
sizes these experiments are run at; dotted segments are hypothetical extensions
that nobody has run — 2,000 inbred lines, 3,000 genotyped outbred individuals.
The eight-way RIL panel is not expected to exceed 750 lines and X-QTL is drawn to
114,000, this experiment. Points mark each design as currently practised, and the
dashed horizontal line is 80% power.

## The point

X-QTL is the least efficient of these designs per fly, but longevity is phenotyped
in bulk, so it is the only one that can be run at 10^5 individuals — and at the
effect sizes resolved here that is the whole of the difference.
