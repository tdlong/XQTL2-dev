#!/bin/bash
# submit.final_rerun.sh — FINAL clean run at the chosen filters, standard pipeline
# end to end (post #22/#23/#24/#25), one command:
#   bash scripts_oneoffs/AGE_SY/catalog_v7/submit.final_rerun.sh
#
# 1. rebuild the catalog at min-dp 10, maxaf 0.03, snpgap 20 (BAQ-on discovery #23,
#    multiallelic positions dropped #24, 1 core #25)  -> process/AGE_SY_v7/Catalog
# 2. count the samples + STANDARD merge (works now #24)  -> process/AGE_SY_v7/Calls/RefAlt.*
# 3. compare_refalt_calls.R v3 vs v7 (SNP counts + overlap + shared-SNP agreement)
#    -> logs/AGE_SY/compare_v3_v7.out
set -euo pipefail
CATDIR=process/AGE_SY_v7/Catalog
DIR=process/AGE_SY_v7
BAMLIST=helpfiles/AGE_SY/bam_list.v4.txt
V3DIR=process/AGE_SY_v3
SELF=scripts_oneoffs/AGE_SY/catalog_v7

git -C pipeline pull --ff-only origin main
echo "pipeline: $(git -C pipeline log --oneline -1)"
for f in "$BAMLIST" pipeline/helpfiles/B_founders.bams.txt; do
  [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

LAST=$(bash pipeline/scripts/build_catalog.sh \
        --founders pipeline/helpfiles/B_founders.bams.txt --out "$CATDIR" \
        --min-dp 10 --maxaf 0.03 --snpgap 20 | tail -1)
echo "catalog build terminal job: $LAST"
[[ "$LAST" =~ ^[0-9]+$ ]] || { echo "ERROR: no job id from build_catalog (got '$LAST')" >&2; exit 1; }

J=$(sbatch --parsable --dependency=afterok:"$LAST" -A tdlong_lab -p standard --time=00:30:00 \
     -o logs/AGE_SY/final_launch.out \
     --wrap="bash $SELF/launch_count_compare.sh $DIR $CATDIR $BAMLIST $V3DIR")
echo "count+compare launcher: $J -> logs/AGE_SY/compare_v3_v7.out (at end)"
echo
squeue -u "$USER"
