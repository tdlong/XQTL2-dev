#!/bin/bash
###############################################################################
# 02_refalt_haps_scan.sh — AGE_SY round 3 (v3): REFALT -> haps -> scan
#
# Runs the full recall/haplotype/scan chain for R1-R12 into process/AGE_SY_v3/.
#   round 1 -> process/AGE_SY      (R1-R7)
#   round 2 -> process/AGE_SY_v2   (R1-R11)
#   round 3 -> process/AGE_SY_v3   (R1-R12)   <-- this script
# A new dir preserves the earlier results.
#
# PREREQUISITES (run 01_align.sh first, then do these before this script):
#   1. R12 BAMs aligned and >1 GB in data/bam/AGE_SY/
#   2. helpfiles/AGE_SY/AGE_SY.bams  includes the 6 R12 sample BAMs
#      (only those passing the 1 GB filter; insert above the founder rows)
#   3. helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R  names_in_bam includes R12
#      (regenerate with the find -size +1G one-liner in README.processing.txt)
#   4. helpfiles/AGE_SY/AGE_SY{10,20}_{F,M}.test.txt  regenerated with R12
#      via ../common/make_AGE_SY_design_files.py (REPS now 1..12; R12 must be
#      present in summary_info_v1.xlsx first)
#
# Upstream barcode / readname mapping (used by 01_align.sh):
#   helpfiles/AGE_SY/readname.mapping.AGE_SY.July_26.txt
#
# Recipe: pipeline/README.md — "run REFALT + haplotypes + multiple scans all at
#         once" (End-to-end section). REFALT + haps run ONCE, then one run_scan.sh
#         per design (4 scans off the same haplotypes). This is the documented
#         multi-scan pattern; run_full_pipeline.sh would redo REFALT+haps per scan.
#         Resources are owned by the run_*.sh wrappers (no #SBATCH here).
#
# Run from: XQTL2-dev repo root, on the cluster.
###############################################################################
set -euo pipefail

DIR=process/AGE_SY_v3
BAMLIST=helpfiles/AGE_SY/AGE_SY.bams
PARFILE=helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R

mkdir -p "$DIR"

JID_REFALT=$(bash pipeline/scripts/run_refalt.sh --bamlist "$BAMLIST" --dir "$DIR")
echo "REFALT: $JID_REFALT"

JID_HAPS=$(bash pipeline/scripts/run_haps.sh --after "$JID_REFALT" --parfile "$PARFILE" --dir "$DIR")
echo "HAPS:   $JID_HAPS"

for scan in AGE_SY10_F AGE_SY10_M AGE_SY20_F AGE_SY20_M; do
    bash pipeline/scripts/run_scan.sh --after "$JID_HAPS" --dir "$DIR" \
        --scan "$scan" --design "helpfiles/AGE_SY/${scan}.test.txt"
    echo "SCAN submitted: $scan"
done

echo ""
echo "All round-3 (v3) jobs submitted -> $DIR"
