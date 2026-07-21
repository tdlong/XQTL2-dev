# AGE_SY_v4 — evaluating the founder-catalog REFALT caller

## Goal
Test the XQTL2 candidate caller (founder-catalog, README appendix "Proposed
founder-catalog REFALT pipeline") on AGE_SY, and compare its callset against the
validated QUAL-based caller (`process/AGE_SY_v3`). Motivation: the prior
investigation ([[NOTES.depth_flip_investigation]] in ../round3_v3_R12) proved the
tiled-vs-whole-chr disagreements (0.02% of sites) are all the pooled `QUAL>59`
filter flipping — never the counts — driven by the `-d` cap and by BAQ moving
QUAL window-dependently. The catalog caller drops QUAL entirely: it ascertains
SNPs once from the founders and counts each BAM at fixed positions with `-B`
(BAQ off), so site membership is deterministic and interval-independent.

## Design (why "first 6 replicates")
AGE_SY = 6 temporal replicates R1-R6 (this batch) + R7-R12 (added later). Each
replicate = 2 sexes (F/M) x 3 treatments (Con, SY10, SY20) = 6 BAMs. Founders =
8 B-population lines (B1-B7, AB8), read from the parfile `founders=` vector.

- **Phase 1**: `bam_list.v4.txt` = 36 samples (R1-R6) + 8 founders = 44 BAMs.
- **Phase 2**: point `--bamlist` at the full `AGE_SY.bams` (80) with the SAME
  `--dir`; catalog + 44 prior counts reused, only 36 new (R7-R12) counted.
  `process/AGE_SY_v4` then has all 72 samples — matches AGE_SY_v3 for compare.

## Files
- `helpfiles/AGE_SY/bam_list.v4.txt` — phase-1 subset (44 BAMs), committed.
- `run_catalog_v4.sh` — orchestrator; `PHASE=1` (default) or `PHASE=2`.
- Config: `helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R` (unchanged; founders + step/size/h_cutoff).

## Run
```bash
bash scripts_oneoffs/AGE_SY/catalog_v4/run_catalog_v4.sh           # phase 1
# after phase 1 completes:
PHASE=2 bash scripts_oneoffs/AGE_SY/catalog_v4/run_catalog_v4.sh   # phase 2
```
`run_refalt.catalog.sh` self-submits its SLURM array + merge and prints the merge
JID. Deliverable: `process/AGE_SY_v4/RefAlt.<chr>.txt`; catalog at
`process/AGE_SY_v4/catalog.tsv.gz` (built on the phase-1 run).

## Compare (the test)
```bash
module load R
Rscript pipeline/scripts/compare_refalt_calls.R process/AGE_SY_v3 process/AGE_SY_v4
```
Pass criteria (README): candidate keeps >= as many usable SNPs with high overlap;
shared-site counts agree closely (residual = BAQ-on v3 vs BAQ-off v4, small mean
|freq diff|); current-only sites are QUAL-passed-but-founder-unclean (few, suspect);
candidate-only sites are real founder-segregating SNPs QUAL dropped (the gain).
NOTE: for a per-sample count comparison, use phase 2 (same 72 samples as v3);
phase 1 differs in sample set, so compare SNP-site membership there, not per-sample counts.

## Status
- [x] Phase-1 bam_list built + verified (36 R1-R6 + 8 founders, no R7-R12 leakage).
- [ ] Phase 1 submitted / done.
- [ ] Phase 2 submitted / done.
- [ ] compare_refalt_calls.R run; result recorded here.

## README-accuracy check (for a possible XQTL2 issue)
Verified the appendix matches the wrapper's actual flags (`--bamlist/--parfile/
--dir/--catalog/--after`) and that the catalog is documented to live in the
project `--dir`, built if absent (README lines ~1108, 1120-1121, 1112-1114) — this
matches how XQTL2 described it. No doc discrepancy found so far; if the run writes
the catalog anywhere other than `--dir`, that's a code/README mismatch -> file an
XQTL2 issue with the evidence.
