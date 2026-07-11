#!/bin/bash
###############################################################################
# run_snpcall_haps.sh — AGE_SY round 3 (v3): fire SNP calling + haplotypes
#
# ONE command to trigger the two multi-day jobs as a dependency chain:
#   REFALT (SNP calling, ~2 days)  ->  haplotypes (~1 day, --after REFALT)
# Both submit immediately and return; haps waits on REFALT via SLURM.
# NEITHER needs the R12 fly counts (N) — those are only for the scans.
#
# Does, in order:
#   1. guard: all six R12 BAMs present and >1 GB (fq2bam finished)
#   2. 03_update_helpfiles.sh: rebuild AGE_SY.bams + names_in_bam from disk
#   3. submit REFALT              -> process/AGE_SY_v3
#   4. submit haps --after REFALT -> process/AGE_SY_v3
#
# Runtime:  REFALT ~2 days, then haps ~1 day. Walk away for the weekend.
#
# Run from: XQTL2-dev repo root, on the cluster (after `git pull`).
###############################################################################
set -euo pipefail
HERE=scripts_oneoffs/AGE_SY/round3_v3_R12
DIR=process/AGE_SY_v3
BAMDIR=data/bam/AGE_SY

# --- 1. guard: R12 alignment must be fully done (protects the 2-day REFALT) ---
missing=0
for sex in F M; do for trt in AgeSY10 AgeSY20 Con; do
    f="${BAMDIR}/${trt}_R12_${sex}.bam"
    if [ -z "$(find "$f" -size +1G 2>/dev/null)" ]; then
        echo "  NOT READY (missing or <1 GB): $f"; missing=1
    fi
done; done
if [ "$missing" != 0 ]; then
    echo "ERROR: not all six R12 BAMs are ready — is fq2bam still running (squeue -u \$USER)?" >&2
    echo "       Aborting so REFALT does not run on incomplete data." >&2
    exit 1
fi

# --- 2. rebuild the bam list + names_in_bam (no fly counts needed) ---
bash "${HERE}/03_update_helpfiles.sh"

# --- 3. SNP calling (REFALT) ---
JID_REFALT=$(bash pipeline/scripts/run_refalt.sh \
    --bamlist helpfiles/AGE_SY/AGE_SY.bams --dir "$DIR")
echo "REFALT (SNP calling) submitted: ${JID_REFALT}   (~2 days)"

# --- 4. haplotypes, chained after REFALT ---
JID_HAPS=$(bash pipeline/scripts/run_haps.sh \
    --after "$JID_REFALT" \
    --parfile helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R --dir "$DIR")
echo "haplotypes submitted:           ${JID_HAPS}   (starts after REFALT, ~1 day)"

cat <<EOF

================================================================
Both jobs are queued. Monitor:  squeue -u \$USER
Output dir:  ${DIR}

Commit the regenerated helpfiles when convenient:
  git add helpfiles/AGE_SY/AGE_SY.bams helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R
  git commit -m "AGE_SY: add R12 to bam list + names_in_bam" && git push

------------------------- NEXT (Monday) ------------------------
The four scans are the ONLY step that needs the R12 fly counts (N).
Once you have them:

  1. Add R12 rows (Num, Percent_Aged) to helpfiles/AGE_SY/summary_info_v1.xlsx
  2. Regenerate the design files:
       python3 scripts_oneoffs/AGE_SY/common/make_AGE_SY_design_files.py
  3. Run the four scans (REFALT+haps will be done by then):
       bash ${HERE}/run_scans.sh

Results land in ${DIR}/<scan>/<scan>.scan.txt
================================================================
EOF
