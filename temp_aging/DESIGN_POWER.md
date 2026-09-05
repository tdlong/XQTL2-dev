# What conventional designs could have detected

4 Sep 2026. Three alternative designs on one common framework, to put a number
behind the claim that regions contributing less than 2% would have been missed
by prior work. Autosomes only.

## Framework

Non-centrality is `lambda = n R^2 / (1 - R^2)`, where `R^2` is the locus
contribution as a fraction of the variance of whatever unit is measured —
individuals in A, line means in B and C. Power is
`P(chi2_df(lambda) > threshold)` at a Bonferroni threshold of 0.05 over the
stated number of tests.

## Base population — outbred diploid, MAF 0.15

| | |
|---|---|
| genotype frequencies | aa 0.0225, Aa 0.2550, AA 0.7225 |
| `Va(locus) = 2pq alpha^2` | 0.02, so `alpha` = 0.280 |
| `Vg` = 0.50 | locus 0.02 + background 0.48 |
| `Ve` | 0.50 |
| `Vp` | 1.00 |

The locus explains 2% of individual phenotypic variance and `h2` is 50%.

## Inbred lines — only aa and AA, at 15% and 85%

| | |
|---|---|
| locus variance among lines | `pq(2 alpha)^2` = 0.0400 — exactly twice the outbred 0.0200 |
| background `Vg` among lines | `2 x 0.48` = 0.96 — unlinked genes double as well |
| total genetic among lines | 1.00, i.e. `2 x Vg`, as it must be |
| `Ve` on a mean of 10 measures | 0.50/10 = 0.050 |
| `Var(line mean)` | 1.050 |
| **locus `R^2`** | **0.0381** |

The causative site does contribute more under inbreeding, but the background
rises by the same factor, so **the doubling cancels in the ratio**. Locus `R^2`
goes from 2.00% in outbred individuals to only 2.67% in unreplicated inbred
lines (0.04/1.5). It reaches 3.81% because replication cuts `Ve` tenfold.
Nearly all of the inbred advantage is the replication, not the homozygosity.

## The three designs

| design | n | df | tests | `R^2` | lambda | threshold | **power** |
|---|---|---|---|---|---|---|---|
| A — outbred diploids, 1 measure each | 1000 | 1 | 1e6 | 0.0200 | 20.4 | 29.7 | **0.175** |
| B — inbred lines, 10 measures each | 200 | 1 | 1e6 | 0.0381 | 7.9 | 29.7 | **0.004** |
| C — 8-way RILs, 10 measures each | 750 | 7 | 1e4 | 0.0381 | 29.7 | 36.9 | **0.462** |

Smallest locus detectable at 80% power: **A 3.96%, B 10.4%, C 2.89%** of
phenotypic variance.

To reach 80% power at a 2% locus each design would need roughly 2000
individuals, 1000 lines, or 1100 RILs.

## Simulation check

2000 replicates per design. A is dosage regression on 1000 individuals; B a
1 df contrast on 200 line means; C a genuine 8-level ANOVA on 750 line means
with founder effects scaled to give an among-line variance of 0.04.

| design | analytic | simulated |
|---|---|---|
| A | 0.175 | 0.171 |
| B | 0.004 | 0.007 |
| C | 0.462 | 0.464 |

Analytic and simulated agree, so the models are set up correctly.

## Conclusion

The 2% figure in the Results is defensible. At a locus explaining 2% of
phenotypic variance no alternative design has better than a 46% chance of
finding it, and the inbred-line panel has essentially none.

Two structural points. B is weakest despite its replication, because a 1 df
test against a million-test threshold with only 200 lines is a poor trade. And
C is the nearest competitor precisely because it shares the 8-way, 7 df
structure and a linkage-scale testing burden — the honest comparison.

## Assumptions that would move the numbers

1. "Explains 2%" is read as 2% of individual **phenotypic** variance in the base
   population. If it means 2% of `Vg`, every `R^2` halves.
2. **MAF 0.15** drops out once the locus is stated as a variance fraction. It
   would matter if an effect size were specified instead.
3. The **inbred doubling** applies to B and C and is doing real work.
4. **1e6 tests for B.** A DGRP-style panel tests roughly 2M common SNPs; that
   change moves the threshold by about 1.4 and barely moves power.
