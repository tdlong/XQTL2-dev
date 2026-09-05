# Proposed changes to RESULTS.md — round 3

Rounds 1 and 2 are applied. One item is open: B4, the 2% claim, which you left
in place pending a real calculation. That calculation is now done —
`temp_aging/design_power.R`, written up in `temp_aging/DESIGN_POWER.md`.

---

## B4 — the 2% now has numbers behind it

Your three designs on one framework, with the models derived rather than
asserted and checked against simulation:

| design | n | df | tests | R² | lambda | **power at a 2% locus** | smallest at 80% |
|---|---|---|---|---|---|---|---|
| A — outbred diploids, 1 measure | 1000 | 1 | 1e6 | 0.0200 | 20.4 | **0.175** | 3.96% |
| B — inbred lines, 10 measures | 200 | 1 | 1e6 | 0.0381 | 7.9 | **0.004** | 10.4% |
| C — 8-way RILs, 10 measures | 750 | 7 | 1e4 | 0.0381 | 29.7 | **0.461** | 2.89% |

Analytic against simulated power: 0.175/0.171, 0.004/0.007, 0.461/0.464.

**Your 2% stands.** No alternative design has better than a 46% chance at that
effect size. To reach 80% they would need about 1940 individuals, 1000 lines,
or 1044 RILs.

The structural point worth having in hand if a reviewer pushes: inbreeding
doubles the locus variance but doubles the unlinked background too, so those
cancel. Locus R² goes 2.00% to only 2.67% from homozygosity alone, and reaches
3.81% because ten measures per line cut Ve tenfold. Almost the whole inbred-line
advantage is replication, not homozygosity.

**Proposed** — replace the placeholder parenthetical currently in paragraph 2:

> Regions contributing less than 2% to heritable variation would have gone
> undetected in prior work employing modestly sized mapping panels consisting of
> several hundred RILs (cites), or association type panels consisting of 100s of
> ILs or 1000s of genotyped individuals. On a common framework a locus
> explaining 2% of phenotypic variance would be found with probability 0.46 by
> 750 eight-way RILs phenotyped ten times each, 0.18 by a thousand outbred
> individuals measured once, and 0.004 by two hundred inbred lines phenotyped
> ten times.

`[ ] approve   [ ] deny   [ ] modify:`

Or keep the current parenthetical placeholder and put the table in the
supplement instead.

`[ ] table to supplement, leave the text as it stands`

---

## Assumptions that would move these numbers

1. **"Explains 2%"** read as 2% of individual *phenotypic* variance in the base
   population. If it means 2% of Vg, every R² halves.
2. **MAF 0.15** drops out once the locus is stated as a variance fraction. It
   would matter only if an effect size were specified instead.
3. The **inbred doubling** applies to B and C and is doing real work.
4. **1e6 tests for B.** A DGRP-style panel tests roughly 2M common SNPs; that
   moves the threshold by about 1.4 and barely moves power.
