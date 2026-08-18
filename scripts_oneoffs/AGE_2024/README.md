# AGE_2024 — the Aug 2024 aging pilot, put on the current pipeline

A preliminary aging experiment: 6 cages, **females only, lab food**, one control
and one selected pool per cage. Sequenced Aug–Oct 2024 and analysed with the
pre-XQTL2 code, so nothing about it was comparable to AGE_SY.

This runs it through the current pipeline at AGE_SY's settings so the two can be
put side by side.

## What runs

```bash
bash scripts_oneoffs/AGE_2024/run_all.sh      # submits everything
```

1. **SNP list** copied from `process/AGE_SY/Catalog`. Same eight founders, so the
   same list applies — and copying rather than rebuilding means both experiments
   are counted at identical sites, which is the point. The eight founder count
   files come across too: same BAMs, same catalog, so identical by construction.
2. **SNPs called per sample** for the 12 BAMs in `data/bam/AGE_2024`.
3. **Haplotypes** at 75 kb windows stepping 5 kb.
4. **Scan** smoothed at 100 kb, design `helpfiles/AGE_Aug13_24/Ageing_Aug13.txt`.
5. **AGE_SY re-scanned on replicates 1–6** — four scans, matching the pilot's
   replicate count. Haplotypes already exist, so this is scan-only.

Steps 2–4 are chained on job IDs; step 5 runs immediately alongside.

Then, once the queue drains:

```bash
Rscript scripts_oneoffs/AGE_2024/gather_scans.R
```

which writes one `process/AGE_2024/AGE_2024_vs_AGE_SY.txt.gz` to scp home —
the pilot, the four matched 6-replicate AGE_SY scans, and the four full
12-replicate ones.

## Why these settings

Every parameter matches `AGE_SY_haplotype_parameters_size75k.R` except the sample
names. Same founders, same window, same step, same smoothing, same SNP set. If
the two scans disagree, it is the experiment and not the analysis.

## What is comparable, and what is not

The pilot is females on lab food, so `AGE_SY10_F` and `AGE_SY20_F` are the
like-for-like comparisons; the male scans are context.

Pool sizes favour the pilot slightly — 482–1177 flies selected per cage against
314–834 for the AGE_SY females — so it is not a weaker experiment per replicate,
just a shorter one. Selection was also more variable (4.0–11.3% against a 5.5%
median in AGE_SY), and the first two cages were weakly selected at 11.3% and 8.7%.

AGE_SY shares one control between diets within a sex; the pilot has its own
control per cage, which is marginally cleaner.

## Files

| file | what |
|---|---|
| `run_all.sh` | submits steps 1–5 |
| `make_reps1_6_designs.R` | cuts the four AGE_SY designs to REP ≤ 6 |
| `gather_scans.R` | collapses all nine scans to one file (run on HPC3) |

Helpfiles live in `helpfiles/AGE_Aug13_24/` beside the original 2024 configs:
`AGE_2024.bams`, `AGE_2024_haplotype_parameters_size75k.R`, and the existing
`Ageing_Aug13.txt` design, which already has the columns the pipeline needs.
