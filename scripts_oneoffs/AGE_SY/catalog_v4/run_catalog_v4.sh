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

JID=$(bash pipeline/scripts/run_refalt.catalog.sh \
        --bamlist "$BAMLIST" \
        --parfile "$PARFILE" \
        --dir     "$DIR")
echo "submitted; merge job id: $JID"
echo
echo "when done, the deliverable is $DIR/RefAlt.<chr>.txt"
echo "compare against the validated QUAL caller (process/AGE_SY_v3):"
echo "  module load R"
echo "  Rscript pipeline/scripts/compare_refalt_calls.R process/AGE_SY_v3 $DIR"
