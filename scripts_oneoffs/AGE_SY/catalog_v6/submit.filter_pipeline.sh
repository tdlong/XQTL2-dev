#!/bin/bash
# submit.filter_pipeline.sh — the WHOLE analysis, one command:
#   bash scripts_oneoffs/AGE_SY/catalog_v6/submit.filter_pipeline.sh
#
# 1. rebuild catalog with BAQ ON (#23), loosest filter (min-dp 10, maxaf 0.07, snpgap 0)
#    -> catalog.annot.tsv.gz + loosest catalog.tsv.gz (superset of all 12 combos)
# 2. SNP-loss table for the 12 combos (indel 0/10/20/50 x maxaf 3/5/7%)  [after build]
# 3. count the 36 samples once at the loosest catalog (BAQ off, no genotype model, #22)
# 4. concordance vs v3 for each of the 12 combos (>2% / >5% ALT-freq diff)  [after count]
#
# Everything chains by SLURM dependency; two logs land: snp_loss_grid.out (early) and
# concordance_grid.out (end).
set -euo pipefail
CATDIR=process/AGE_SY_v6/Catalog
DIR=process/AGE_SY_v6
BAMLIST=helpfiles/AGE_SY/bam_list.v4.txt
V3DIR=process/AGE_SY_v3
SELF=scripts_oneoffs/AGE_SY/catalog_v6

git -C pipeline pull --ff-only origin main
echo "pipeline: $(git -C pipeline log --oneline -1)"
for f in "$BAMLIST" pipeline/helpfiles/B_founders.bams.txt; do
  [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

# 1. build (BAQ-on) + loosest filter; build_catalog prints its terminal job id
LAST=$(bash pipeline/scripts/build_catalog.sh \
        --founders pipeline/helpfiles/B_founders.bams.txt --out "$CATDIR" \
        --min-dp 10 --maxaf 0.07 --snpgap 0 | tail -1)
echo "catalog build chain terminal job: $LAST"
[[ "$LAST" =~ ^[0-9]+$ ]] || { echo "ERROR: did not get a job id from build_catalog (got '$LAST')" >&2; exit 1; }

# 2. SNP-loss table (after build)
J_LOSS=$(sbatch --parsable --dependency=afterok:"$LAST" -o logs/AGE_SY/snp_loss_grid.out \
        "$SELF/snp_loss_grid.sh" "$CATDIR")
echo "SNP-loss table  : $J_LOSS -> logs/AGE_SY/snp_loss_grid.out"

# 3+4. count + concordance (gate job runs after build, self-submits the rest)
J_GATE=$(sbatch --parsable --dependency=afterok:"$LAST" -A tdlong_lab -p standard \
        --time=20:00 -o logs/AGE_SY/phase2_launch.out \
        --wrap="bash $SELF/launch_count_concordance.sh $DIR $CATDIR $BAMLIST $V3DIR")
echo "count+concordance launcher: $J_GATE -> logs/AGE_SY/concordance_grid.out (at end)"
echo
squeue -u "$USER"
