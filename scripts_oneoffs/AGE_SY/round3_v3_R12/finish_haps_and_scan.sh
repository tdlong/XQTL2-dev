#!/bin/bash
###############################################################################
# finish_haps_and_scan.sh — one submit to run the rest of AGE_SY v3 to completion.
#
# chr3R haplotypes timed out (run_haps.sh overrode the step's 1-day budget with
# 4h). This re-runs ONLY chr3R by calling REFALT2haps.sh directly (so its own
# 1-day header applies), then chains the four scans to start once it finishes.
#
#   chr3R haps (REFALT2haps --array=5, highmem, 1 day)
#      └─ afterok ─> smooth -> R2 -> hap_scan -> concat   x4 designs (run_scans.sh)
#
# The other 4 chromosomes' haplotypes are already done; the scans depend on the
# chr3R haps job so they start only when all 5 are complete. afterok means a
# failed haps job stops the scans instead of scanning incomplete data.
#
# Run once from repo root on the cluster:
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/finish_haps_and_scan.sh
###############################################################################
set -euo pipefail
HERE=scripts_oneoffs/AGE_SY/round3_v3_R12
DIR=process/AGE_SY_v3
PARFILE=helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R

# --- guards: everything the chain needs must already be in place ---
[ -f "$DIR/RefAlt.chr3R.txt" ] || { echo "ERROR: $DIR/RefAlt.chr3R.txt missing (SNP calling incomplete?)" >&2; exit 1; }
for c in chrX chr2L chr2R chr3L; do
  [ -f "$DIR/R.haps.$c.out.rds" ] || { echo "ERROR: haps not done for $c — expected $DIR/R.haps.$c.out.rds" >&2; exit 1; }
done
for s in AGE_SY10_F AGE_SY10_M AGE_SY20_F AGE_SY20_M; do
  [ -f "helpfiles/AGE_SY/$s.test.txt" ] || { echo "ERROR: design missing: helpfiles/AGE_SY/$s.test.txt" >&2; exit 1; }
done
echo "Prereqs OK: RefAlt 5/5, haps 4/5 (chr3R pending), 4 design files present."

# --- 1. finish chr3R haplotypes (direct call = 1-day header, not run_haps.sh's 4h) ---
JID_HAPS=$(sbatch --parsable --array=5 pipeline/scripts/REFALT2haps.sh \
    --parfile "$PARFILE" --dir "$DIR")
echo "chr3R haplotypes submitted: $JID_HAPS  (highmem, 10G, 1 day)"

# --- 2. the four scans, chained to start after chr3R haps completes ---
AFTER="$JID_HAPS" bash "$HERE/run_scans.sh"

cat <<EOF

================================================================
Queued as one dependency chain — nothing else to submit:
  chr3R haps ($JID_HAPS)
     -> smooth -> R2 -> hap_scan -> concat   for each of the 4 designs

Monitor:  squeue -u \$USER
Final results (4 files):  $DIR/<scan>/<scan>.scan.txt
================================================================
EOF
