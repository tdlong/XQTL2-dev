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
- 2026-07-24 eve: SUBMITTED, running overnight.
  - baq_caller (v4 counting-recipe diagnosis): array `54412588_[1-4]` + merge `54412589`.
    -> logs/AGE_SY/baq_caller_merge.out ; process/AGE_SY_v4/baq_caller.Con_R5_F.chrX.tsv.gz
  - v5 catalog build: array `54412606_[1-5]` (per-chr) -> gather `54412607` -> filter
    `54412608` -> cleanup `54412609`. -> process/AGE_SY_v5/Catalog/ (catalog.annot.tsv.gz,
    catalog.tsv.gz, catalog.stats.txt).
  NEXT (when back): `bash scripts/cluster_sync.sh`, then read baq_caller_merge.out
  (BAQ vs -C alleles verdict) and catalog.stats.txt (min-dp10/snpgap25 SNP counts vs v4).
  Do NOT resubmit either chain — both are already queued/running.

## 2026-07-27 RESULTS

### baq_caller (v4 counting-recipe diagnosis) — recipe is NOT the cause
All 4 tasks agree (array ran redundantly; merge crashed on a --basename wrinkle,
irrelevant — each task printed the full summary):
  catalog(-B,-C)=32.070  +BAQ(-C)=31.549  +free(-B,-m)=32.069  v3ish(-m)=31.549  (ALT%)
- -C alleles constraint: ZERO effect (32.070 vs 32.069).
- -q/-Q: zero (already shown).
- BAQ: the ONLY lever, and only 0.5pp (~1.6%), BAQ-ON (=v3) *removes* alt.
Direction+size are OPPOSITE and ~4x smaller than the ~6.5% "v4 drops alt" gap.
=> Re-counting the SAME bam at the SAME positions v3-flags vs v4-flags gives ~identical
counts. The counting recipe does NOT explain the v3-vs-v4 disagreement.
=> Prime suspect is the CONFIRMED #14 dup-position bug in v4's RefAlt (delta_gt2 was
built on that table). #14 is fixed in the pipeline v5 is built on. DECISIVE TEST =
call v5 clean, redo v3-vs-v5 count compare; if the gap vanishes, it was #14.

### v5 catalog build — clean
candidate 2,000,582 | near-indel(snpgap25) -426,704 | depth(min-dp10) -187,052 |
not-near-fixed -134,337 | not-segregating -67,553 | KEPT 1,184,936.
min-dp10 no longer the throttle; snpgap=25 is now the biggest single cut (427k) ->
MEASURE it from the dist-to-indel annotation (catalog.annot.tsv.gz, re-cut in seconds).

### NEXT: v5 call step
`bash scripts_oneoffs/AGE_SY/catalog_v5/submit.call_v5.sh` — 36 samples (same v4 bam
list) against v5 catalog; then v3-vs-v5 count compare.

## 2026-07-27 THE ANSWER — v4 loses ALT because `call -m -C alleles` genotypes it away

Dissection path (one sample, one chr, real reads — not statistics):
- Con_R5_F chrX, filtered mid-arm 5-18Mb + cov 50-300x + gap>10% = **3,537** disagreeing SNPs.
- Founder triage: 493 (14%) near-indel (filter them, separate issue); **3,044 (86%) clean**,
  of which **2,961 are v4-loses-ALT** — the dominant phenomenon.
- Dissected 5 clean v4-loses-ALT SNPs (dissect_snp.sh flag ladder). Decisive:
  the ALT reads survive ALL quality filters (mpileup `-B -q20 -Q20` keeps 33/26/37/35 ALT),
  then `| bcftools call -m -C alleles` sets ALT="." and AD collapses to REF-only (0 ALT).
  Only diff between the two = the `call` step. Example chrX:6471332 A,G: mpileup 212,33 ->
  post-call 227,0. (5th SNP 5332918 ~18% ALT cleared the caller threshold, kept.)

MECHANISM: per-sample **diploid genotype model** in `call -m` calls hom-REF at a
low-but-real pooled ALT fraction (~13%) and deletes the allele. v3 escapes it by calling
JOINTLY across all 36 samples (cohort-confident). This is the misapplied prior TL flagged
from the start: at a KNOWN pooled SNP we want to COUNT ref vs alt, not genotype.

FIX (filed): count AD straight from `bcftools mpileup -a FORMAT/AD`, never route through
`call -m -C alleles`. Draft issue in-repo: `issue_catalog_count_call.md` (file with gh).
This is THE cause of the v4/catalog undercount at clean SNPs, distinct from the near-indel
(BAQ, v3-side) tail.
