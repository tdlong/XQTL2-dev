# AGE_SY catalog v5 — lab notebook

## Why v5
Pipeline advanced to `5843d33` (XQTL2 #21: annotate SNPs, threshold downstream;
#20: good_SNPs segregation-test fix). v4 was built on the older caller. v5 rebuilds
on the new one for the phase-2 comparison.

**What the new caller changes:** *which SNPs* are in the catalog — build now emits
`catalog.annot.tsv.gz` (every candidate SNP + distance-to-indel + per-founder AD),
and `catalog_filter.sh` cuts thresholds from it in seconds (no rebuild). New defaults:
`--min-dp 10` (was 20, #15), `--snpgap 25` (was 5, #21), `--maxaf 0.03`,
`--exempt-founders B5:chr2L`.

**What it does NOT change:** *how reads are counted* per SNP. `catalog_count.sh`'s
`-q20 -Q20 -B --max-depth 2000` recipe is untouched — that's the separate layer the
`baq_caller` experiment (v4 dir) is diagnosing (BAQ vs `-C alleles` as the source of
v4's ALT loss vs v3). The two are orthogonal and run in parallel.

## Plan
1. **Build** (this step, ~1.5 h): `submit.build_v5.sh` → `process/AGE_SY_v5/Catalog/`.
   Checkpoint: eyeball `catalog.stats.txt` (per-rule SNP tally under min-dp 10 /
   snpgap 25) and compare kept-SNP counts to v4 before spending the call step.
2. **Call + verify** (next): 36 samples, same v4 BAM list
   (`helpfiles/AGE_SY/bam_list.v4.txt`) so v4-vs-v5 differ ONLY in caller version.
3. Later: re-cut snpgap from the dist-to-indel annotation (measure the right width);
   phase-2 add R7–R12 → 72 samples.

## Log
- 2026-07-24: confirmed cluster pipeline at 5843d33, catalog_filter.sh present.
  Built `submit.build_v5.sh` (pull --ff-only + build_catalog with new defaults).
  baq_caller array (54412588) + merge (54412589) running in parallel on v4.
