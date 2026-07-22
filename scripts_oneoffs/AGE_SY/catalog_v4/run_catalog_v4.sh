#!/usr/bin/env bash
# run_catalog_v4.sh — evaluate the proposed founder-catalog REFALT caller on AGE_SY.
#
# Runs the XQTL2 candidate caller (pipeline/scripts/run_refalt.catalog.sh, README
# appendix "Proposed founder-catalog REFALT pipeline") into process/AGE_SY_v4,
# then we compare it against the validated QUAL-based calls in process/AGE_SY_v3.
#
# The candidate defines the SNP set ONCE from the 8 B-founders (built into the
# --dir on first run) and counts every BAM against that fixed catalog with -B
# (BAQ off) at fixed -T positions — deterministic, interval-independent, and
# incremental (adding a sample = one more count job, nothing recalled).
#
# Two phases (the point of the test is the incremental machinery):
#   PHASE 1 (this script, default): 6 temporal replicates R1-R6
#           = 3 treatments x 2 sexes x 6 = 36 samples + 8 founders (44 BAMs).
#           Builds catalog, counts 44 BAMs, merges RefAlt.<chr>.txt.
#   PHASE 2 (PHASE=2): add replicates R7-R12 -> all 72 samples. Point --bamlist
#           at the full AGE_SY.bams (80 BAMs) with the SAME --dir; the catalog and
#           the 44 already-counted BAMs are reused, only the 36 new are counted.
#           process/AGE_SY_v4 then matches AGE_SY_v3's 72 samples for a clean compare.
#
# run_refalt.catalog.sh self-submits its SLURM array + merge and prints the merge
# JID, so this is a login-node orchestrator (bash it, do NOT sbatch it).
#
# Run from repo root:
#   bash scripts_oneoffs/AGE_SY/catalog_v4/run_catalog_v4.sh          # phase 1
#   PHASE=2 bash scripts_oneoffs/AGE_SY/catalog_v4/run_catalog_v4.sh  # phase 2 (after phase 1 done)
set -uo pipefail

PARFILE=helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R   # founders=c("B1"..."AB8")
DIR=process/AGE_SY_v4
PHASE=${PHASE:-1}

if [[ "$PHASE" == "1" ]]; then
  BAMLIST=helpfiles/AGE_SY/bam_list.v4.txt               # 36 samples (R1-R6) + 8 founders
  echo "PHASE 1: catalog build + count 44 BAMs (R1-R6 + founders) -> $DIR"
elif [[ "$PHASE" == "2" ]]; then
  BAMLIST=helpfiles/AGE_SY/AGE_SY.bams                   # all 72 samples + 8 founders
  echo "PHASE 2: incremental — add R7-R12 (36 new counts), reuse catalog + 44 prior -> $DIR"
else
  echo "unknown PHASE='$PHASE' (use 1 or 2)" >&2; exit 1
fi

for f in "$BAMLIST" "$PARFILE" pipeline/scripts/run_refalt.catalog.sh; do
  [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done

# ---------------------------------------------------------------------------
# PRE-FLIGHT: say out loud what will be REUSED vs COMPUTED, so the final
# RefAlt.<chr>.txt's provenance is never a mystery. The pipeline silently reuses
# any output file that already exists in --dir; this makes that visible up front.
# Run with DRYRUN=1 to print this and stop WITHOUT submitting anything.
# ---------------------------------------------------------------------------
echo "== pre-flight: $DIR =="
npieces=$(ls "$DIR"/catalog.chr*.bed 2>/dev/null | wc -l | tr -d ' ')
if [[ "${npieces:-0}" -ge 1 ]]; then
  echo "  catalog build   : REUSE $npieces existing catalog.chr*.bed piece(s) — founders NOT recalled"
else
  echo "  catalog build   : COMPUTE from founders (SLOW, ~1.5 h)"
fi
if [[ -s "$DIR/catalog.tsv.gz" ]]; then
  echo "  catalog assemble: REUSE existing catalog.tsv.gz  (first line: $(zcat "$DIR/catalog.tsv.gz" 2>/dev/null | head -1))"
else
  echo "  catalog assemble: COMPUTE (assemble from pieces)"
fi
nbam=$(grep -cve '^[[:space:]]*$' "$BAMLIST")
ncount=$(ls "$DIR"/counts/*.tsv.gz 2>/dev/null | wc -l | tr -d ' ')
echo "  sample counts   : REUSE ${ncount:-0} already-counted, COMPUTE $((nbam - ${ncount:-0})) new (of $nbam BAMs)"
echo "  -> to force any of these fresh, DELETE the file(s) above first; rerunning alone never recomputes them."
echo
if [[ "${DRYRUN:-0}" == "1" ]]; then
  echo "DRYRUN=1 set — not submitting. Rerun without DRYRUN to launch."
  exit 0
fi

JID=$(bash pipeline/scripts/run_refalt.catalog.sh \
        --bamlist "$BAMLIST" \
        --parfile "$PARFILE" \
        --dir     "$DIR")
echo "submitted; merge job id: $JID"

# Chain the verification (dup check + per-chr counts vs v3 + compare_refalt_calls.R
# + B5 depth) to run automatically after the merge succeeds, so this is one kickoff.
SELFDIR=$(cd "$(dirname "$0")" && pwd)
VJID=$(sbatch --parsable --dependency=afterok:"$JID" "$SELFDIR/compare_and_diagnose.sh")
echo "verification chained after merge: job $VJID -> logs/AGE_SY/compare_v3_v4_${VJID}.out"
echo
echo "deliverable: $DIR/RefAlt.<chr>.txt"
echo "when everything finishes:  bash scripts/cluster_sync.sh"
