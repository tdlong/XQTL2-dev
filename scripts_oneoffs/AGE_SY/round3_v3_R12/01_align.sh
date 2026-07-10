#!/bin/bash
###############################################################################
# 01_align.sh — AGE_SY round 3 (v3): align R12 reads (July_26 sequencing run)
#
# Round 3 adds replicate R12 to the aging experiment. R1-R11 BAMs already
# exist in data/bam/AGE_SY/. R12 was a SINGLE sequencing run, so there is
# NO merge step here (contrast round 2 / R8-R11, which came in two runs and
# needed ../round2_v2_R8-R11/mergebams.sh).
#
# ---------------------------------------------------------------------------
# BARCODE / READNAME MAPPING (the demux table: <i7>  <i5>  <sample>)
#     helpfiles/AGE_SY/readname.mapping.AGE_SY.July_26.txt
#   6 rows -> {Con,AgeSY10,AgeSY20}_R12_{F,M}
#   (all barcode maps for this project live in helpfiles/AGE_SY/ as
#    readname.mapping.AGE_SY.<RUN>.txt — one file per sequencing run.)
# ---------------------------------------------------------------------------
#
# RAW READS:
#   /dfs7/adl/tdlong/fly_pool/XQTL2-dev/data/raw/AGE_SY/July_26
#   filenames: xR096-L1-G5-P0nn-<i7>-<i5>-R1/R2.fastq.gz
#
# OUTPUT BAMs -> data/bam/AGE_SY/   (final dir; no temp/merge needed)
#
# Recipe: pipeline/README.md — "Adding replicates to an existing experiment",
#         Step 1 (align new samples only). fq2bam invocation per §"Step 2 —
#         Align reads". Resources are owned by fq2bam.sh (no #SBATCH here).
#
# Run from: XQTL2-dev repo root, on the cluster.
###############################################################################
set -euo pipefail

# --- where the barcodes are stored (see banner above) ---
BARCODES=helpfiles/AGE_SY/readname.mapping.AGE_SY.July_26.txt
RAWDIR=/dfs7/adl/tdlong/fly_pool/XQTL2-dev/data/raw/AGE_SY/July_26
BAMDIR=data/bam/AGE_SY

[ -f "$BARCODES" ] || { echo "ERROR: barcode mapping not found: $BARCODES" >&2; exit 1; }
[ -d "$RAWDIR" ]   || { echo "ERROR: raw read dir not found: $RAWDIR" >&2; exit 1; }

N=$(grep -cve '^[[:space:]]*$' "$BARCODES")   # number of samples (6)
echo "Barcode map : $BARCODES  (${N} samples)"
echo "Raw reads   : $RAWDIR"
echo "Output BAMs : $BAMDIR"
echo ""

JOBID=$(sbatch --array=1-${N} --parsable \
    pipeline/scripts/fq2bam.sh "$BARCODES" "$RAWDIR" "$BAMDIR")
echo "fq2bam array job submitted: $JOBID"
echo ""
echo "When it finishes, verify R12 BAM sizes (>1 GB) before 02_refalt_haps_scan.sh:"
echo "  ls -l $BAMDIR/{Con,AgeSY10,AgeSY20}_R12_{F,M}.bam"
