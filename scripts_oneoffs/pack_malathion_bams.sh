#!/bin/bash
# pack_malathion_bams.sh
# Run this on the server to bundle the malathion BAMs for public distribution.
# Output: malathion_bams.tar in the current directory.
#
# Usage (from anywhere):
#   bash /dfs7/adl/tdlong/fly_pool/XQTL2-dev/scripts_oneoffs/pack_malathion_bams.sh

BAMDIR=/dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/data/bam

echo "Packing malathion BAMs from ${BAMDIR}..."
tar -cvf malathion_bams.tar \
    -C "${BAMDIR}" \
    c_F_1.bam c_F_1.bam.bai \
    c_M_1.bam c_M_1.bam.bai \
    s_F_1.bam s_F_1.bam.bai \
    s_M_1.bam s_M_1.bam.bai

echo ""
echo "Done: malathion_bams.tar ($(du -sh malathion_bams.tar | cut -f1))"
echo "Copy to wfitch: scp malathion_bams.tar tdlong@wfitch.bio.uci.edu:~/public_html/"
