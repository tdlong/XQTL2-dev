#!/bin/bash
# run_ZINC2_v4.sh — Rerun ZINC2 haplotype calling + scan with fixed pipeline
#
# Starts from existing REFALT files in process/ZINC2/.
# Runs REFALT2haps (fixed haplotype estimator) then full scan for F and M.
#
# Run from: XQTL2-dev root

set -euo pipefail


# Submit haplotype job once; chain both scans off it
JID=$(bash pipeline/scripts/run_haps.sh \
    --parfile helpfiles/ZINC2/ZINC2_haplotype_parameters.R \
    --dir     process/ZINC2)

bash pipeline/scripts/run_scan.sh \
    --after  ${JID} \
    --dir    process/ZINC2 \
    --scan   ZINC2_F_v4 \
    --design helpfiles/ZINC2/Zinc2.test.F.txt

bash pipeline/scripts/run_scan.sh \
    --after  ${JID} \
    --dir    process/ZINC2 \
    --scan   ZINC2_M_v4 \
    --design helpfiles/ZINC2/Zinc2.test.M.txt
