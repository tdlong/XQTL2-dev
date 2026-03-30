#!/bin/bash
# malathion_pipeline.sh
# Training run: malathion resistance XQTL in D. melanogaster (A-population)
#
# BAMs (c_F_1, c_M_1, s_F_1, s_M_1) must be in data/bam/malathion/
# and indexed (.bam.bai).  Download from:
#   https://wfitch.bio.uci.edu/~tdlong/malathion_bams.tar
# then:
#   mkdir -p data/bam/malathion
#   tar -xf malathion_bams.tar -C data/bam/malathion/
#
# Run from the XQTL2-dev repo root.

set -euo pipefail

PROJECT=malathion
DESIGN=helpfiles/${PROJECT}/design.txt
PARFILE=helpfiles/${PROJECT}/hap_params.R
PROCDIR=process/${PROJECT}
SCAN=malathion_scan

mkdir -p ${PROCDIR}

# ── Step 3: REFALT counts ────────────────────────────────────────────────────
# Build bam list: pooled samples + A-pop founders
find data/bam/${PROJECT} -name "*.bam" -size +100M \
    | sort > helpfiles/${PROJECT}/bam_list.txt
cat pipeline/helpfiles/founder.bams.txt >> helpfiles/${PROJECT}/bam_list.txt

JID_REFALT=$(sbatch --parsable --array=1-5 \
    -A tdlong_lab -p standard \
    pipeline/scripts/bam2bcf2REFALT.sh \
        helpfiles/${PROJECT}/bam_list.txt \
        ${PROCDIR})
echo "Step 3 submitted: ${JID_REFALT}"

# ── Step 4: Call haplotypes ──────────────────────────────────────────────────
JID_HAPS=$(sbatch --parsable --array=1-5 \
    --dependency=afterok:${JID_REFALT} \
    -A tdlong_lab -p highmem --mem-per-cpu=10G \
    pipeline/scripts/REFALT2haps.sh \
        --parfile ${PARFILE} \
        --dir     ${PROCDIR})
echo "Step 4 submitted: ${JID_HAPS}"

# ── Step 5: Scan ─────────────────────────────────────────────────────────────
bash pipeline/scripts/run_scan.sh \
    --design  ${DESIGN} \
    --dir     ${PROCDIR} \
    --scan    ${SCAN} \
    --after   ${JID_HAPS} \
    -A tdlong_lab -p standard

echo "Pipeline submitted. Monitor with: squeue -u \$USER"
