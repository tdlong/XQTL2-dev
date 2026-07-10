#!/bin/bash
# pack_malathion_bams.sh
# Run this on the server to bundle the malathion BAMs for public distribution.
# Output: /dfs7/adl/tdlong/fly_pool/XQTL2/data/malathion_bams.tar
#
# Usage:
#   bash /dfs7/adl/tdlong/fly_pool/XQTL2-dev/scripts_oneoffs/pack_malathion_bams.sh

BAMDIR=/dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/data/bam
OUTFILE=/dfs7/adl/tdlong/fly_pool/XQTL2/data/malathion_bams.tar

mkdir -p "$(dirname ${OUTFILE})"

echo "Packing malathion BAMs from ${BAMDIR}..."
tar -cvf "${OUTFILE}" \
    -C "${BAMDIR}" \
    c_F_1.bam c_F_1.bam.bai \
    c_M_1.bam c_M_1.bam.bai \
    s_F_1.bam s_F_1.bam.bai \
    s_M_1.bam s_M_1.bam.bai

echo ""
echo "Done: ${OUTFILE} ($(du -sh ${OUTFILE} | cut -f1))"
echo "Copy to wfitch: scp ${OUTFILE} tdlong@wfitch.bio.uci.edu:~/public_html/"
