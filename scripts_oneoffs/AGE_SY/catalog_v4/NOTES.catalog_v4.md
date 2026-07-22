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
- [x] Phase 1 — ran after #13 fix (catalog built, 44 BAMs counted, RefAlt produced).
- [x] compare_and_diagnose.sh run (log: logs/AGE_SY/compare_v3_v4_54375203.out).
- [ ] Phase 2 — hold until XQTL2 #14/#15 addressed (no point widening a malformed table).

## RESULT + two XQTL2 issues filed (2026-07-22)
Evaluation of the founder-catalog caller vs v3 surfaced two problems, both filed:

1. **BUG #14 (confirmed, high priority): RefAlt duplicates nearly every position.**
   `RefAlt.chrX.txt` = 196,427 lines but only 66,810 unique POS; `uniq -d` on
   (CHROM,POS) = **66,793 duplicated positions**. Root cause: `catalog_count.sh:67`
   writes each sample's OWN observed REF/ALT (not the catalog's fixed REF/ALT), and
   `catalog_merge.R:33` outer-joins on (CHROM,POS,REF,ALT) -> positions fragment
   when ALT text varies per sample. Corrupts the count table + the v3-vs-v4 compare
   (the "0% identical / 0.35 freq diff" is this artifact) and would break REFALT2haps.
   Fix: emit/merge on the catalog's fixed REF/ALT (or merge on CHROM,POS only).
   https://github.com/tdlong/XQTL2/issues/14

2. **DESIGN #15: `--min-dp 20` (all founders) too aggressive -> lower to ~10x.**
   v4 keeps ~1/4 the SNPs of v3 (chrX 66,810 vs 250,801; chr2L 15,208 vs 326,237).
   Per-founder mean DP (chr2R): B1 175, B3 121, AB8 119, B7 61, B2 36, B4 34,
   B6 21, B5 15 (B5 chr2L 7.3). Requiring all 8 >= 20x at every site is the throttle;
   chr2L gutted by B5. B5 is fixed-but-shallow (the `*.RG.bam`, reference reads
   subtracted to fix an inversion het), so a lower floor is safe.
   https://github.com/tdlong/XQTL2/issues/15

Both belong to XQTL2 — do NOT patch pipeline scripts here (CLAUDE.md rule 2).
NOTE: earlier "compare on quality" pass-criteria numbers are invalid until #14 is
fixed (the comparison table is corrupted by the duplication).

## BLOCKER: XQTL2 issue #13 (module htslib clash)
Phase-1 catalog_build failed immediately (exit 127, 0-byte founders.calls.bcf).
Root cause: `catalog_build.sh` (and `catalog_count.sh`) load BOTH `bcftools/1.21`
and `samtools/1.10`; samtools/1.10's htslib/1.10.2 shadows the htslib>=1.16 that
bcftools/1.21 needs -> every bcftools call dies with
`undefined symbol: bcf_has_variant_types, version HTSLIB_1.16`. This is the SAME
clash documented in ../round3_v3_R12/NOTES.depth_flip_investigation.md.
Filed: https://github.com/tdlong/XQTL2/issues/13 (with error, root cause, and 3
fix options: read SM tag via `bcftools view -h` / isolate samtools in a subshell /
use samtools built on htslib>=1.16). NOT our bug — do not patch pipeline scripts
here (CLAUDE.md rule 2). Resume phase 1 once XQTL2 ships the fix and we pull it.

## README-accuracy check (for a possible XQTL2 issue)
Verified the appendix matches the wrapper's actual flags (`--bamlist/--parfile/
--dir/--catalog/--after`) and that the catalog is documented to live in the
project `--dir`, built if absent (README lines ~1108, 1120-1121, 1112-1114) — this
matches how XQTL2 described it. No doc discrepancy found so far; if the run writes
the catalog anywhere other than `--dir`, that's a code/README mismatch -> file an
XQTL2 issue with the evidence.
