#!/bin/bash
# run_ZINC2_v4.sh — Rerun ZINC2 haplotype calling + scan with fixed pipeline
#
# Starts from existing REFALT files in process/ZINC2/.
# Runs REFALT2haps (fixed haplotype estimator) then full scan for F and M.
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2

set -euo pipefail

PARFILE=helpfiles/ZINC2/ZINC2_haplotype_parameters.R
SNP_TABLE=pipeline/helpfiles/FREQ_SNPs_Bpop.cM.txt.gz
FOUNDERS=B1,B2,B3,B4,B5,B6,B7,AB8

# Run haplotypes + female scan; capture haplotype job ID
F_OUT=$(bash pipeline/scripts/run_full_pipeline.sh \
    --project      ZINC2 \
    --parfile      ${PARFILE} \
    --design       helpfiles/ZINC2/Zinc2.test.F.txt \
    --scan         ZINC2_F_v4 \
    --snp-table    ${SNP_TABLE} \
    --founder-list ${FOUNDERS} \
    --skip-fq2bam --skip-refalt)
echo "$F_OUT"
JID_HAPS=$(echo "$F_OUT" | grep "^haplotypes:" | awk '{print $2}')

# Male scan reuses same haplotypes — no recomputation
bash pipeline/scripts/run_full_pipeline.sh \
    --project      ZINC2 \
    --parfile      ${PARFILE} \
    --design       helpfiles/ZINC2/Zinc2.test.M.txt \
    --scan         ZINC2_M_v4 \
    --snp-table    ${SNP_TABLE} \
    --founder-list ${FOUNDERS} \
    --skip-fq2bam --skip-refalt --skip-haps --after ${JID_HAPS}
