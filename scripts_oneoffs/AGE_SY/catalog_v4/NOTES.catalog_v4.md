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

## UPDATE 2026-07-22: XQTL2 fixed #13, #14, #15 (all CLOSED)
- #14 fix (commit 851f37b): catalog_count.sh now `bcftools call -m -C alleles -T cat`
  (constrains every sample to the catalog's REF/ALT; emits CHROM POS + AD only);
  catalog_merge.R keys on (CHROM,POS) alone. -> one row per position, no fragmentation.
- #15 fix: MIN_DP default 10 (was 20), PLUS new `--exempt-founders` (default
  `B5:chr2L`) dropping B5 from the depth/fixation gate on chr2L only. Rationale:
  B5 chr2L is non-polymorphic BY CONSTRUCTION (reads forced to an ALT-only ref in
  the inversion fix), so founder rules don't apply there. Added data/founders/FOUNDERS.md.
- Our run_catalog_v4.sh unchanged (exemption default baked into catalog_build.sh).

RERUN TO CONFIRM (stale process/AGE_SY_v4 must be wiped — catalog+count logic changed):
  rm -rf process/AGE_SY_v4
  bash scripts_oneoffs/AGE_SY/catalog_v4/run_catalog_v4.sh
  bash scripts/cluster_sync.sh
Check: (1) 0 duplicate positions, (2) chr2L recovered toward ~326k, (3) site counts up.
- [x] rerun attempted (job 54377xxx).

## RERUN 1 FAILED — #14 fix incomplete (REOPENED #14) 2026-07-22
Build (54377096_1-5) + gather COMPLETED, but ALL 44 catalog_count tasks FAILED
(exit 127, <1s); merge never ran (DependencyNeverSatisfied). Log:
  Unable to parse the -T file; expected CHROM\tPOS\tREF,ALT with -C alleles
  but found instead: chr2L 9204 T C
Cause: the #14 fix changed count to `bcftools call -m -C alleles -T "$cat"`
(catalog_count.sh:68), which needs the catalog as 3-col CHROM POS REF,ALT
(REF,ALT comma-joined). But catalog_build.sh:146 still writes 4-col CHROM POS
REF ALT. Fix updated the consumer, not the catalog format.
Already FIXED upstream by 227200a (catalog_gather.sh now emits CHROM POS REF,ALT
comma-joined) — landed 15:39, ~11 min AFTER our count failed (15:28). Our run used
the pre-fix pipeline + a stale 4-col catalog.tsv.gz. (I reopened #14 on stale info;
corrected + re-closed — the format is in catalog_gather.sh, not catalog_build.sh.)

RESUME WITHOUT REBUILD (the slow catalog.<chr>.bed pieces are good and reused;
only the seconds-long gather produced the wrong-format tsv):
  git -C pipeline pull origin main
  rm -f process/AGE_SY_v4/catalog.tsv.gz process/AGE_SY_v4/catalog.tsv.gz.tbi
  rm -rf process/AGE_SY_v4/counts
  bash scripts_oneoffs/AGE_SY/catalog_v4/run_catalog_v4.sh   # build tasks skip; gather regenerates; count+merge
  bash scripts/cluster_sync.sh
- [ ] duplication gone / counts recovered (rerun with fixed gather).

## RERUN 2 — catalog correct, but two more issues filed 2026-07-22
Catalog regenerated to correct comma format (chr2L 9204 T,C). Kicked count+merge;
user cancelled after seeing a tabix error in a count log. Two issues filed:
- #16: catalog_count.sh prints `[E::get_intv] Failed to parse ... "CHROM POS
  REF_x ALT_x"` on EVERY count, yet the job succeeds (unskipped header in
  `tabix -s1 -b2 -e2`). Result is fine (merge reads whole file, .tbi unused) but a
  clean run emitting an ERROR makes logs untrustworthy. Fix: `tabix -S 1`.
  https://github.com/tdlong/XQTL2/issues/16
- #17: catalog steps silently reuse any existing output by FILENAME (no provenance,
  no staleness check) -> a final RefAlt can silently mix stale intermediates and its
  origin is unknowable. Fix: loud reuse logging + provenance manifest + hash-based
  staleness + explicit --fresh. https://github.com/tdlong/XQTL2/issues/17
Our run_catalog_v4.sh now prints an explicit REUSE-vs-COMPUTE pre-flight (+DRYRUN)
and auto-chains the verification after merge.
- #18: catalog_count.sh over-requests 2 cores + 12GB; seff shows ~1 core (52% of 2)
  and 268MB used. One count task per sample (44+/run) -> wasted core-hours. Fix:
  --cpus-per-task=1. https://github.com/tdlong/XQTL2/issues/18
  (seff also CONFIRMS the count job COMPLETED exit 0 despite the #16 tabix error.)
Open XQTL2 issues from this evaluation: #16 (tabix log), #17 (provenance), #18 (cores).

OPEN (my offer): make run_catalog_v4.sh run COUNT-ONLY when the catalog exists
(catalog_count.sh array + merge), so the build never re-queues. Awaiting user call.

## REWORKED WORKFLOW 2026-07-24 (XQTL2 main @ d40e279) — two explicit commands
XQTL2 split the caller into build vs call, with an organized output tree. Fixes
along the way: #14 (comma catalog), #16 (dropped vestigial count tabix), #18 (count
1 core). run_refalt.catalog.sh RETIRED. New layout:
  process/AGE_SY_v4/Catalog/  catalog.tsv.gz (+ catalog.stats.txt per-rule tally)
  process/AGE_SY_v4/Calls/    RefAlt.<chr>.txt (+ counts/)

Two commands (build once, then call/add samples):
  # 1. build catalog from B-pop founders (ONE-TIME, ~1.5 h, OVERWRITES --out)
  bash pipeline/scripts/build_catalog.sh \
      --founders pipeline/helpfiles/B_founders.bams.txt --out process/AGE_SY_v4/Catalog
  # 2. call samples (+ founders as columns) and auto-verify (our wrapper)
  bash scripts_oneoffs/AGE_SY/catalog_v4/run_catalog_v4.sh
  # 3. push logs
  bash scripts/cluster_sync.sh

- B_founders.bams.txt == our 8 founders incl B5.RG.bam (matches bam_list.v4.txt).
- call_samples counts EXACTLY --bamlist (no skip); founders in the list => columns
  (matches v3). Add samples later: rerun run_catalog_v4.sh with BAMLIST=just-new-BAMs.
- Our scripts updated for the new layout: run_catalog_v4.sh now wraps call_samples
  (+errors if catalog not built) and chains compare_and_diagnose; compare_and_diagnose
  reads Calls/ + Catalog/catalog.stats.txt.
- User already rm -rf process/AGE_SY_v4 (fresh start).
- [x] catalog built + checked (catalog_check_54411533.out).

## CATALOG CHECK (54411533) — HEALTHY, chr2L recovered
Rules: 1,775,608 candidates -> 1,324,782 KEPT (75%). Drops: founder-depth 218k
(12%), not-near-fixed 158k (9%), not-segregating 74k (4%). Format: chr2L 9204 T,C
(comma-joined, #14 baked in). Per-chr catalog vs v3 RefAlt:
  chrX 176,670/250,802 (0.70); chr2L 298,812/326,238 (0.92 — RECOVERED from 15k);
  chr2R 244,953/348,347 (0.70); chr3L 304,342/411,873 (0.74); chr3R 300,005/406,693 (0.74).
chr2L fixed by B5:chr2L exemption + min-dp 10. Catalog keeps 70-92% of v3 sites
(FEWER than v3, vs README "keeps >= as many") — expected for founder-segregating
vs QUAL>59; resolve whether v3-only sites are junk-or-real via compare's v3-only
output AFTER calling. VERDICT: catalog good -> proceed to call R1-R6.
- [ ] call R1-R6 + verify (dup 0, counts vs v3, v3-only sites junk?).

## KEY INSIGHT (reported to dev #1): catalog == v3's downstream good_SNPs filter
The catalog's ascertainment rules ARE the `good_SNPs` filter the current pipeline
already applies DOWNSTREAM in REFALT2haps.code.R:169-175:
  zeros==0        = every founder covered      = catalog min-dp
  notfixed==0 @0.03/0.97 = near-fixed per founder = catalog maxaf 0.03
  informative     = segregating                = catalog polymorphic/segregating
So v3 RefAlt = RAW (QUAL>59); the SNPs actually USED = v3-AFTER-good_SNPs. The catalog
IS that filtered set, ascertained early (no QUAL). => v4 keeping 70-92% of RAW v3 is
EXPECTED; the fair comparison is v4 vs v3-AFTER-good_SNPs.
PROPOSED TEST (next analysis): apply good_SNPs to v3's RefAlt (has founder cols),
compare to v4 catalog. Match => catalog reproduces the pipeline's founder-clean set
(the goal). Differ => that delta is exactly what QUAL was adding/removing.
Reported: https://github.com/tdlong/XQTL2-dev/issues/1#issuecomment-5074550560

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
