# AGE_SY — founder-catalog caller evaluation

Master record for the evaluation of XQTL2's founder-catalog SNP caller against the
validated QUAL caller, on AGE_SY. One place; the per-phase `NOTES.*.md` hold detail.

## Goal
Decide whether the founder-catalog caller is trustworthy for pooled REF/ALT counting,
and pick the catalog filters. Baseline for comparison = the validated QUAL caller.

## Directories (process/)
| dir | what it is | status |
|-----|------------|--------|
| `AGE_SY_v3` | validated QUAL caller output (`RefAlt.*`, `calls.*.bcf`, `R.haps`) | **BASELINE — keep** |
| `AGE_SY_v6` | BAQ-on catalog (annot) + 44-sample counts + `snp_pass` + `merge_dedup` RefAlt | **LIVE — keep** |
| ~~`v4`~~ | first eval (buggy #14 counts) + `baq_caller`/`flag_ladder`/`delta_gt2` diagnostics | deleted |
| ~~`v5`~~ | BAQ-off catalog build, never counted | deleted |
| ~~`v7`~~ | wasteful final-rerun stub (cancelled) | deleted |

## Script folders (this dir; `common`/`getdata`/`round*` are the experiment proper, not this eval)
| folder | role | status |
|--------|------|--------|
| `catalog_v6/` | filter-analysis pipeline (snp_loss + concordance) + `merge_dedup` **hack** | live |
| `catalog_v7/` | subset-existing-counts + compare_refalt_calls | live |
| ~~`catalog_v4`~~ | the diagnosis (agreement_delta, flag_ladder, baq_caller, dissect → #22) | deleted; in git history |
| ~~`catalog_v5`~~ | build tooling, founder_triage, filter_grid, dissect_snp | deleted; in git history |

When #26 lands and the clean run replaces the `merge_dedup` hack, `catalog_v6` + `catalog_v7`
collapse into a single `catalog/` folder.

## What we found (→ XQTL2 issues)
- Counts were being corrupted / deflated, traced to specific causes, each filed:
  - **#14** RefAlt duplicate positions (fixed) · **#15** min-dp too aggressive · **#16–19** misc.
  - **#20** good_SNPs no-op · **#21** annotate-once / filter-downstream design.
  - **#22** `catalog_count` ran `call -m -C alleles`, whose diploid genotype model **deletes
    real minor-allele reads** in pooled samples — the core bug. Fixed: count AD from mpileup.
  - **#23** BAQ was off during founder discovery (should be on); fixed.
  - **#24** `catalog_merge` cartesian-explodes on multiallelic positions (hard blocker); fixed
    at build. **#26** = move that drop downstream so it doesn't force a founder recall + document
    rerun impact. **#25** catalog_build CPU over-allocation.

## Decision (see `catalog_v6/NOTES.catalog_v6.md` for the 12-combo table)
**min-dp 10, maxaf 3%, snpgap 20.** Indel distance is the lever (knee at 20); maxaf is a
non-lever; BAQ-on discovery already removed most near-indel junk.

## Current state / next step
Waiting on **#26** (multiallelic drop as a downstream filter). Once it lands, do the clean
final run **in `v6`, in place** (no new `vN`): `catalog_filter` re-cut of the existing annot
(no founder recall) → standard `catalog_count` + `catalog_merge` → `compare_refalt_calls.R`
vs v3. Then the `merge_dedup` hack is deleted and `catalog_v6`/`catalog_v7` collapse into one.
