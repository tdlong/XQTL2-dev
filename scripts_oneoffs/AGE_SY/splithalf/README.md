# AGE_SY split-half — odd vs even replicates

**Not the main pipeline.** Nothing here writes inside `process/AGE_SY`. The live
scans and figures cannot be overwritten by anything in this folder.

## Why

To partition the spatial variation in Cutler H² into a mean, a sex effect, a
sugar effect, and their interaction, we need to know how much of what we see is
estimation noise. Running each of the four treatments twice — once on odd
replicates, once on even — gives that directly: noise is independent between the
two halves, so across windows the **covariance** between the odd and the even
estimate of any contrast is signal, and the shortfall against its total variance
is noise. No modelling assumptions, and it works even though the four contrasts
do not carry equal noise (SY10 and SY20 share their control samples, so control
noise cancels in the sugar and interaction contrasts but not in the sex one).

Odd/even rather than R1–R6 vs R7–R12 because the replicates arrived over several
months — a block split would put batch straight into the error term.

## What does not survive this

SY20 selected a larger fraction of flies than SY10 (males 0.068 vs 0.054 mean,
females 0.055 vs 0.045). That is systematic, so it sits in both halves and reads
as signal, not noise. Only matters if a sugar main effect turns up.

## Layout

```
helpfiles/AGE_SY/splithalf/         8 design files, .odd / .even
process/AGE_SY_splithalf/           everything this produces
  Haps -> ../AGE_SY/Haps            symlink; haplotypes reused, not recomputed
  Scans/AGE_SY10_F_odd/ …           8 scans
  AGE_SY_splithalf_H2.txt.gz        the long file
```

## Running it

From the repo root on the cluster:

```bash
Rscript scripts_oneoffs/AGE_SY/splithalf/make_splithalf_designs.R
bash    scripts_oneoffs/AGE_SY/splithalf/run_splithalf_scans.sh
```

The second submits all 8 scans and chains the gather job behind them, so the
long file appears on its own. Log: `logs/AGE_SY/splithalf_gather.<jobid>.out`.

`smooth=100` kb, matching the live run (haplotypes are `size=75000`).

## The output

`process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz`, tab-separated:

| column | |
|---|---|
| `chr`, `pos` | window, `pos` in bp |
| `sex` | `F`, `M` |
| `sugar` | `SY10`, `SY20` |
| `half` | `odd`, `even` |
| `Cutl_H2` | Cutler broad-sense heritability |
| `Falc_H2`, `Wald_log10p` | carried along; free, and save a second trip |

One row per window per cell, 8 cells. Small enough to pull down and explore
locally, which is the point:

```bash
scp hpc3:<path>/process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz .
```

## Analysis (local, once the long file is down)

```bash
Rscript scripts_oneoffs/AGE_SY/splithalf/partition_H2.R        # the numbers
Rscript scripts_oneoffs/AGE_SY/splithalf/plot_H2_components.R  # the picture
```

`partition_H2.R` reports, per term, the signal variance, the noise variance, the
odd-even reliability with a block-bootstrap interval, and the share of signal.

`plot_H2_components.R` splits the heritability at each position into its mean,
sex, diet and sex:diet components: the total signal along the genome, and
beneath it the share attributable to each term. Two things make it honest —

- each component is `contrast_odd * contrast_even`, not the square. A square is
  noise-inflated, so at a background window all four squares are noise of
  similar size and the stack shows a tidy 25/25/25/25 that looks like a result.
  The product's noise cancels in expectation, so background sits at zero.
- shares are drawn only where the total clears `median + 3 MAD`. A quantile cut
  admits noise, and noise over noise fills the panel with vivid nonsense.

Both were checked on synthetic data with a sex effect planted on chrX and an
interaction on chr3L: the components peak at 0.0395 and 0.0532 at the planted
positions, against 0.0005-0.0007 for the two terms with nothing in them.

## Depends on XQTL2 #32

The design files keep their real `REP` values — the odd file is `1,3,5,7,9,11`.
Before #32, `Heritability()` matched replicates by position rather than by
label, so those files would have silently produced H² from 3 of the 6
replicates with mismatched selection proportions, while `Wald_log10p` looked
perfectly healthy. Do not run this against a pipeline predating that fix.
