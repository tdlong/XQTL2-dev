# AGE_SY — founder-catalog caller evaluation

Master record for the evaluation of XQTL2's founder-catalog SNP caller against the
validated QUAL caller, on AGE_SY. Detail in `catalog_v6/NOTES.md`.

## Data lineage (process/ folders) — history, then the final two
- `AGE_Aug13_24` — early **pilot, different raw data** (not an AGE_SY version; stands alone). **KEPT.**
- `AGE_SY` (50) → `AGE_SY_v2` (74) → `AGE_SY_v3` (80 = 72 samples + 8 founders) — the *same*
  experiment re-analyzed as replicates arrived (R1–R6 … → R7–R12). `v3` was the complete
  QUAL-caller baseline; the earlier two superseded partials. **All removed.**
- `AGE_SY_v3_size75k` — v3 at a 75 kb window (symlinked counts). **Removed.**
- `AGE_SY_v6_size75k` — the **new catalog caller** on the full dataset, validated to reproduce v3.
  **Renamed to `AGE_SY`** — the go-forward result.

**Final state: two folders — `AGE_Aug13_24` (pilot) and `AGE_SY` (the new catalog caller).**
The removed QUAL data is regenerable from BAMs, and its validation record lives in git
(`logs/AGE_SY/compare_v3_v6.out`, `logs/figures/`, `catalog_v6/NOTES.md`). Note: figure files
and `catalog_v6/NOTES.md` still say `AGE_SY_v6_size75k` — that folder is now `AGE_SY`.

## Goal
Decide whether the founder-catalog caller is trustworthy for pooled REF/ALT counting,
and pick the catalog filters. Baseline = the validated QUAL caller.

## Result (validated)
The founder-catalog caller reproduces the validated caller to **~0.2% median frequency
difference** at usable depth (≥100×), **no rare-allele bias** (the old genotype-collapse
bug is gone), the only systematic being a ≤0.35% BAQ effect at common alleles — and it
additionally **recovers ~25K chr2L SNPs** the QUAL caller missed. Filters:
**min-dp 10, maxaf 3%, snpgap 20** (indel distance is the lever, knee at 20; maxaf a
non-lever). See `catalog_v6/NOTES.md` for the 12-combo table and the comparison.

**End-to-end (plug-and-play):** called reps 7–12, then ran the existing haplotype+scan
pipeline UNCHANGED on the v6 RefAlt — the size75k genome scan
(`AGE_SY_v6_size75k_4scan.png`) reproduces the July-15 v3 scan essentially exactly (same
QTL peaks/heights, incl. the chr3L ~200 peak), with marginally more chrX signal. The caller
is a drop-in from counts through QTL. → Recommend it become the XQTL2 default.

## Directories (process/)
| dir | what it is | status |
|-----|------------|--------|
| `AGE_Aug13_24` | early pilot, different raw data | keep (distinct experiment) |
| `AGE_SY` | the **new catalog caller** — full R1–R12 catalog + counts + 75 kb haps + scans (was `AGE_SY_v6_size75k`) | **LIVE** |
| ~~`AGE_SY`(old), `_v2`, `_v3`, `_v3_size75k`, `_v4/5/6/7`~~ | superseded QUAL partials + eval iterations | deleted (regenerable; validation in git) |

## Script folder — `catalog_v6/` (single, live)
(`common`/`getdata`/`round*` are the experiment proper, not this eval.)
- `submit.filter_pipeline.sh` + `snp_loss_grid.sh` + `concordance_grid.*` +
  `launch_count_concordance.sh` — build the BAQ-on catalog and sweep filters (the decision).
- `submit.clean_v6.sh` + `launch_count_compare.sh` — the clean final call: `catalog_filter`
  re-cut of the annot → standard count + merge → compare vs v3. No hacks.
- `compare_summary.R` + `submit.compare_summary.sh` — overlap counts + |freq diff| vs
  coverage/MAF (ECDF plots to `logs/figures/`).
- `NOTES.md` — the lab notebook.
- Earlier `catalog_v4`/`v5`/`v7` folders and the `merge_dedup`/`subset` hacks were deleted
  once the pipeline fixes made the standard path work; they live in git history.

## What we found (→ XQTL2 issues, all fixed)
- **#14** RefAlt duplicate positions · **#15** min-dp too aggressive · **#16–19** misc ·
  **#20** good_SNPs no-op · **#21** annotate-once / filter-downstream.
- **#22** `catalog_count` ran `call -m -C alleles` — its diploid genotype model **deleted
  real minor-allele reads** in pooled samples (the core bug). Fixed: count AD from mpileup.
- **#23** BAQ off during founder discovery (should be on) · **#25** catalog_build CPU.
- **#24** `catalog_merge` cartesian-explodes on multiallelic catalog positions ·
  **#26** moved that drop downstream (no founder recall to re-cut) ·
  **#27** `catalog_count` missing `-I` → SNP+indel double-rows crashed the merge.

## State
**Done.** Investigation complete, caller validated, tree tidied. `v3` (baseline) + `v6`
(live) are the only AGE_SY eval dirs; `catalog_v6/` is the only script folder.
