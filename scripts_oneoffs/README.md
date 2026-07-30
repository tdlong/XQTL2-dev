# scripts_oneoffs/ — historical record, not runnable tooling

**Assume nothing in here runs as-is.** These are the submission scripts as they were
when each analysis was done. They are kept so we can see exactly what was run, not so
they can be re-run.

Two pipeline changes have since moved out from under them:

1. **The repo split** (pipeline moved to XQTL2, reached via the `pipeline/` symlink).
   Older scripts call `scripts/run_scan.sh`, `scripts/haps2scan.freqsmooth.sh`,
   `scripts/concat_Chromosome_Scans.Andreas.sh`, `scripts/plot_pseudoscan.R` — bare
   `scripts/` paths from when pipeline and project were one repo. Several of those
   files no longer exist in any form.
2. **XQTL2 #28** — the founder-catalog caller became the default SNP caller, and
   pipeline output moved into stage subfolders (`Catalog/ Calls/ Haps/ Scans/`).
   The joint QUAL caller moved to `pipeline/scripts/legacy/`, so call sites for
   `bam2bcf2REFALT.sh`, `run_refalt.sh`, `run_refalt.tiled.sh`, and
   `compare_refalt.sh` point at paths that no longer resolve. Scripts that build
   their own internal paths (`SCAN_DIR=${DIR}/${SCAN}`, preflight checks on
   `$DIR/R.haps.*.out.rds`) are looking in the pre-#28 flat locations.

These were deliberately **not** patched. Repointing a script we have no intention of
running buys nothing and makes the record less faithful to what was actually executed.

## If you want to re-run one of these analyses

Port it, don't patch it: write a fresh script against the current pipeline —
`build_catalog.sh` + `call_samples.sh` for Step 3, and the `Catalog/Calls/Haps/Scans`
layout throughout. Read `pipeline/README.md` first. The old script is useful as a
statement of the design and parameters, not as a starting point for edits.

Before porting, check whether you need to re-run at all — the results may still be on
the cluster:

```bash
bash scripts/process_fingerprint.sh <prefix>   # e.g. AGE_SY, ZINC, PUPATION, MALATHION
bash scripts/cluster_sync.sh                   # bring the fingerprint back
```

A project reported as `layout: flat (pre-#28)` has usable results but the current
pipeline cannot read them; migrate it in place (no recompute) with
`pipeline/scripts/reorganize_project.sh`. **A project reported as `layout: MIXED` must
not be handed to that script as-is** — it does a bare `mv RefAlt.*.txt Calls/`, so
top-level copies (often stale symlinks) overwrite the real merged tables in `Calls/`.
Remove the redundant top-level copies first.

## Current state by project

- **AGE_SY** — already on the new caller. `process/` holds `AGE_Aug13_24` (early pilot,
  different raw data) and `AGE_SY` (the catalog caller on the full R1–R12 dataset,
  formerly `AGE_SY_v6_size75k`). The superseded QUAL dirs were deleted as regenerable;
  the validation record is in git. See `AGE_SY/README.md`. Nothing to port.
- **ZINC2, ZINC_Hanson, PUPATION, MALATHION** — not audited since #28. Fingerprint them
  before assuming anything about what still exists or which layout it is in.
- **`_shared/`** — dead against the current pipeline on both counts above.
