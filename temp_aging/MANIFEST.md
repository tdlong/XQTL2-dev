# temp_aging — what is here and where it came from

AGE_SY lifespan XQTL: figures, the numbers in them, and the scripts that produce
both. Written so that no file in this folder has to be reverse-engineered.

**Read this first.** Every number currently in the prose comes from the
**12-replicate** dataset (replicates 1–12). The current dataset is the
**10-replicate** one (replicates 8 and 9 dropped — they are the May 2023 cage,
see `helpfiles/AGE_2024/population_assignment.txt`). Only `make_figure1.R` and
`make_figure2.R` can produce it; every other script and every prose file is still
12-replicate. Status is marked per file below.

## The data these read

Nothing here computes haplotype frequencies or scans; all of that is on HPC3.
Five files are fetched by scp and everything below reads them.

| file | made by | what it is |
|---|---|---|
| `process/AGE_SY/AGE_SY_4scan.txt.gz` | `make_4scan_df.R` on HPC3 | the four 12-replicate scans, one row per window × treatment |
| `process/AGE_SY/AGE_SY_4scan_no89.txt.gz` | `scripts_oneoffs/AGE_SY/nov_only/gather.R` | the same at 10 replicates |
| `process/AGE_SY_splithalf/AGE_SY_splithalf_H2.txt.gz` | `scripts_oneoffs/AGE_SY/splithalf/gather_splithalf_H2.R` | eight scans, 4 treatments × odd/even halves, 6+6 |
| `process/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz` | `nov_only/gather.R` | the same at 5+5 |
| `process/AGE_SY/AGE_SY_zoom_means.txt.gz` | `make_zoom_means.R` on HPC3 | per-replicate founder frequencies, ±0.6 Mb around seven peaks |

Derived locally:

| file | made by |
|---|---|
| `process/AGE_SY_splithalf/H2_varcomp_by_window.txt.gz` | `scripts_oneoffs/AGE_SY/splithalf/varcomp_H2.R` |
| `…_no89.txt.gz` | the same, with both paths passed as arguments |

## Figures

| script | reads | writes | status |
|---|---|---|---|
| `make_figure1.R` | 4scan, splithalf H2, varcomp | `Figure1_plot.png`; `Figure1b_plot.png` with `no89` | **current** — takes a variant argument |
| `make_figure2.R` | 4scan, zoom means | `Figure2_plot.png`; `Figure2b_plot.png` with `no89` | **current** — filters replicates 8, 9 locally |

Both are the only files here that handle the 10-replicate data.

## Cluster extracts

| script | runs on | writes | status |
|---|---|---|---|
| `make_4scan_df.R` | HPC3 | `AGE_SY_4scan.txt.gz` | 12-replicate only; the 10-replicate equivalent is `nov_only/gather.R` |
| `make_zoom_means.R` | HPC3 | `AGE_SY_zoom_means.txt.gz` | **current** — extracts all replicates, so both variants use it |

## Analysis — each backs a specific claim

All five print to stdout only; nothing is saved, so a number in the prose is
reproduced by re-running the script. **All five are hardcoded to the
12-replicate files.**

| script | reads | backs |
|---|---|---|
| `h2_threshold.R` | 4scan | the h²-at-Wald-5 table in `results_draft.md` — "windows above ~0.8% of phenotypic variance are significant" |
| `scan_resolution.R` | 4scan | paragraph 1 of `results_draft.md` — the twelve peaks above −log10 P 15 and their 2-unit support intervals |
| `significant_regions.R` | 4scan, splithalf H2 | paragraph 3 — the partition over significant tiles, and the male-vs-female medians by arm |
| `partition_by_wald_rank.R` | splithalf H2 | **half-live.** Its *ratio* column (sex+diet as a share, 9.4% [5.9, 15.9]) is quoted in `results_draft.md` and is sound — numerator and denominator run over the same windows, so the fifteenfold window overlap cancels. Its *absolute* column, rescaling the genome to total 50%, is **withdrawn**: summing overlapping windows and rescaling to the quantity being estimated is circular. `results_draft.md` now states plainly that no genome-wide total is claimed. `SUMMARY.md` still carries the old table and contradicts it. |

## Prose

| file | status |
|---|---|
| `METHODS.md` | mostly replicate-count agnostic — design, the Cutler estimator, the split-half error term, blocking. Needs the replicate counts changed and the rescaling paragraph cut. |
| `results_draft.md` | three paragraphs plus the numbers behind them and a section on what is deliberately not claimed. All 12-replicate. |
| `Figure1_legend.md` | Tony's wording. Numbers are 12-replicate. |
| `Figure2_legend.md` | Tony's wording. Carries a flag that the seven loci are not the seven strongest — how they were chosen still needs stating. |
| `SUMMARY.md` | the long write-up. Contains the withdrawn 50%-rescaled table. **Most stale of the six.** |
| `FLOOR_PROBLEM.md` | about the pipeline's `Cutl_H2_bias` term, not about replicate count. Current, and worth reporting upstream to XQTL2. |

## archive/ — untracked

`pilot_AGE_2024/` the August 2024 pilot and the May/Nov cage follow-up.
`superseded/` scripts replaced by better ones. Both have their own README.

## What has to happen for the 10-replicate dataset

1. Give the five analysis scripts the same variant argument the figures have.
2. Re-run them against `_no89` and renumber the prose. Known shifts: autosomal
   sex 7.3 → 7.8%, diet 2.2 → 3.9%, chr3L peak 207.5 → 163.8.
3. Cut the rescaled column from `partition_by_wald_rank.R` and from `SUMMARY.md`.
4. Say in `Figure2_legend.md` how the seven loci were chosen.
