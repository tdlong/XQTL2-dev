#!/bin/bash
###############################################################################
# run_scans.sh — AGE_SY round 3 (v3): the four scans (Monday step)
#
# Run AFTER run_snpcall_haps.sh has finished (REFALT + haplotypes complete in
# process/AGE_SY_v3) AND the design files include R12's fly counts.
#
# PREREQ (the one thing that needs the R12 N's):
#   1. add R12 rows to helpfiles/AGE_SY/summary_info_v1.xlsx
#   2. python3 scripts_oneoffs/AGE_SY/common/make_AGE_SY_design_files.py
#
# Submits SY10/SY20 x F/M, each smoothing -> Wald+H2 -> concat (pipeline
# run_scan.sh). Haplotypes already exist, so no --after is needed; set
# AFTER=<jobid> in the environment if you want to chain onto a running haps job.
#
# Runtime: hours.  Run from: XQTL2-dev repo root, on the cluster.
###############################################################################
set -euo pipefail
DIR=process/AGE_SY_v3
AFTER_FLAG=""
[ -n "${AFTER:-}" ] && AFTER_FLAG="--after ${AFTER}"

for scan in AGE_SY10_F AGE_SY10_M AGE_SY20_F AGE_SY20_M; do
    design="helpfiles/AGE_SY/${scan}.test.txt"
    [ -f "$design" ] || { echo "ERROR: design not found: $design (run make_AGE_SY_design_files.py)" >&2; exit 1; }
    bash pipeline/scripts/run_scan.sh $AFTER_FLAG \
        --dir "$DIR" --scan "$scan" --design "$design"
    echo "SCAN submitted: $scan"
done
echo ""
echo "All four scans submitted -> ${DIR}/<scan>/<scan>.scan.txt"
