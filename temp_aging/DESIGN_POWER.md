# What conventional designs could have detected

Supporting the claim in Results paragraph 2 that regions contributing less than
2% of variation would have gone undetected in prior work. Autosomes only.

Computed by `temp_aging/design_power.R`; the figure is
`temp_aging/make_FigS2.R` -> `temp_aging/FigureS2_plot.png`.

## The point

Per fly, X-QTL is much **less efficient** than any design that genotypes and
phenotypes individuals. Pooling and sequencing only the selected tail discards
what a replicated line mean keeps: 750 RILs phenotyped ten times each is 7,500
flies for power 0.46, while X-QTL screening the same 7,500 flies gets 0.008.

But the flies are nearly free. The trait is phenotyped **in bulk** — no
individual is ever measured, one lets the cage age and collects the survivors —
and genotyping is on pools rather than individuals.

What limits the other designs is not statistical efficiency but mechanics. Lines
have to be constructed and maintained; individuals have to be handled, genotyped
and phenotyped one at a time. Those panels are more efficient per fly and are
struggling to scale for exactly that reason — a RIL panel is not expected to
exceed 750 lines however large the budget. Bulk phenotyping and pooled
sequencing remove that barrier, which is why 10^5 individuals is routine here
where 10^3 is a major undertaking elsewhere.

Note the asymmetry in the figure. The X-QTL point at 114,000 flies per treatment
is solid: it is an experiment that has been run. Every extension that would let
another design approach it — 2,000 inbred lines, 3,000 genotyped outbred
individuals — is dotted, and there is no particular reason to expect those curves
to be realised.

**So the design is inefficient per fly and efficient per unit of effort, and it
is the second that decides what can be detected.** The condition is that the
trait admits bulk selection. Longevity does: survivors select themselves. For a
trait that requires measuring individuals one at a time the argument collapses
and the conventional designs are the better choice.

## Framework

Non-centrality is `lambda = n R^2 / (1 - R^2)`, where `R^2` is the locus
contribution as a fraction of the variance of whatever unit is measured —
individuals in A, line means in B and C. Power is
`P(chi2_df(lambda) > threshold)` at each design's own genome-wide bar.

Base population: `Vp = 1`, `h2 = 50%` (so `Vg = Ve = 0.5`), MAF 0.15 for the two
SNP designs. Inbred lines and RILs are phenotyped ten times each.

### Base population — outbred diploid, MAF 0.15

| | |
|---|---|
| genotype frequencies | aa 0.0225, Aa 0.2550, AA 0.7225 |
| `Va(locus) = 2pq alpha^2` | 0.02, so `alpha` = 0.280 |
| `Vg` = 0.50 | locus 0.02 + background 0.48 |
| `Ve` | 0.50 |

### Inbred lines — only aa and AA, at 15% and 85%

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
lines (0.04/1.5); it reaches 3.81% because replication cuts `Ve` tenfold.
Nearly all of the inbred-line advantage is the replication, not the homozygosity.

## The four designs

Thresholds: 0.05/1e6 for the SNP panels, 0.05/1e4 for the RIL linkage scan, and
the paper's own LWP 7.5 for X-QTL — stricter than a Bonferroni over 268 tiles,
so X-QTL is handicapped here rather than flattered.

| design | n | df | tests | power at 2% | power at 0.5% |
|---|---|---|---|---|---|
| A — outbred panel, 1 measure | 1,000 | 1 | 1e6 | 0.175 | 0.0007 |
| B — inbred lines, 10 measures | 200 | 1 | 1e6 | 0.004 | 0.0000 |
| C — 8-way RILs, 10 measures | 750 | 7 | 1e4 | 0.461 | 0.0039 |
| D — X-QTL, this experiment | 114,000 | 7 | — | 1.000 | 0.894 |

Smallest locus detectable at 80% power, at the size each design is run at:
**A 3.96%, B 10.4%, C 2.89%**, against **0.42%** for X-QTL at its LWP 7.5
threshold. To reach 80% power at a 2% locus the conventional designs would need
roughly 1,940 individuals, 1,000 lines, or 1,044 RILs.

## The figure

`FigureS2_plot.png` plots power against individuals/lines phenotyped, two
panels: a 2% locus and a 0.5% locus, the size this experiment actually resolves.

Curves stop where each design realistically stops. Solid is the size these
experiments are run at; dotted is a plausible extension:

- **inbred lines** 200 → 2,000, the rumoured extension of the panel
- **outbred panel** 1,000 → 3,000, a guess, nothing supports a specific ceiling

Both are labelled *hypothetical* in the figure legend. Neither has been run, and
the figure should not be read as predicting that either will be.
- **8-way RILs** no extension, the panel is not expected to exceed 750
- **X-QTL** no extension, 114,000 is this experiment per treatment

The left panel is the honest counterweight — at a 2% locus the RIL panel does
respectably and X-QTL needs 30,000 flies to match 750 RILs. The right panel is
the argument: at 0.5% nothing else lifts off the floor even along its dotted
extension (inbred lines reach ~0.15 at 2,000, the outbred panel ~0.05 at 3,000)
while X-QTL is at 0.89.

## These are idealised power calculations for the conventional designs

Each of A, B and C assumes a single homogeneous panel measured under one
protocol: no batch structure, no heterogeneity between sub-experiments, every
individual contributing full information to one analysis. Large phenotyping
efforts are usually assembled over several experiments, and a panel of *n*
assembled that way delivers less power than a panel of *n* run at once.

The X-QTL curve is not idealised in the same way. Its non-centrality comes from
the fitted relation `h2 = 0.00774 (T - 7)`, estimated from the observed scan, so
replicate structure, smoothing and haplotype reconstruction error are already
inside it.

**The comparison therefore favours the conventional designs**, and the gap in
the figure is if anything conservative. This is a general property of textbook
power calculations and needs no claim about any particular published study.

## Simulation check

2000 replicates per design, at the 2% locus. A is dosage regression on 1,000
individuals; B a 1 df contrast on 200 line means; C a genuine 8-level ANOVA on
750 line means with founder effects scaled to an among-line variance of 0.04.

| design | analytic | simulated |
|---|---|---|
| A | 0.175 | 0.171 |
| B | 0.004 | 0.007 |
| C | 0.461 | 0.464 |

Analytic and simulated agree, so the models are set up correctly.

## Assumptions that would move the numbers

1. "Explains 2%" is read as 2% of individual **phenotypic** variance in the base
   population. If it means 2% of `Vg`, every `R^2` halves.
2. **MAF 0.15** drops out once the locus is stated as a variance fraction. It
   would matter if an effect size were specified instead.
3. The **inbred doubling** applies to B and C and is doing real work.
4. **1e6 tests for B.** A DGRP-style panel tests roughly 2M common SNPs; that
   moves the threshold by about 1.4 and barely moves power.
5. **X-QTL held to LWP 7.5.** A Bonferroni over 268 tiles would be 1.9e-4, far
   laxer, and would put its smallest detectable locus at 0.24% rather than 0.42%.
