#!/bin/bash
# submit.clean_v6.sh — the CLEAN final run, in place in v6, standard pipeline only
# (post #22/#23/#24/#25/#26). One command:
#   bash scripts_oneoffs/AGE_SY/catalog_v6/submit.clean_v6.sh
#
# 1. catalog_filter on the EXISTING v6 annot: drop multiallelic (#26) + min-dp10/maxaf3/
#    snpgap20 -> clean catalog.tsv.gz. Seconds, NO founder recall.
# 2. standard catalog_count + catalog_merge -> process/AGE_SY_v6/Calls/RefAlt.* (no hack).
# 3. compare_refalt_calls.R v3 vs v6 -> logs/AGE_SY/compare_v3_v6.out
set -euo pipefail
CATDIR=process/AGE_SY_v6/Catalog
DIR=process/AGE_SY_v6
BAMLIST=helpfiles/AGE_SY/bam_list.v4.txt
V3DIR=process/AGE_SY_v3
SELF=scripts_oneoffs/AGE_SY/catalog_v6

git -C pipeline pull --ff-only origin main
echo "pipeline: $(git -C pipeline log --oneline -1)"
[[ -s "$CATDIR/catalog.annot.tsv.gz" ]] || { echo "no annot at $CATDIR/catalog.annot.tsv.gz" >&2; exit 1; }
[[ -s "$BAMLIST" ]] || { echo "no bamlist $BAMLIST" >&2; exit 1; }

# 1. downstream re-cut on the existing annot (no recall)
JF=$(sbatch --parsable -o logs/AGE_SY/clean_v6_filter.out \
      pipeline/scripts/catalog_filter.sh --catdir "$CATDIR" --min-dp 10 --maxaf 0.03 --snpgap 20)
echo "catalog_filter (re-cut, no recall): $JF -> logs/AGE_SY/clean_v6_filter.out"

# 2+3. after the re-cut: standard count+merge, then compare
JG=$(sbatch --parsable --dependency=afterok:"$JF" -A tdlong_lab -p standard --time=00:30:00 \
      -o logs/AGE_SY/clean_v6_launch.out \
      --wrap="bash $SELF/launch_count_compare.sh $DIR $CATDIR $BAMLIST $V3DIR")
echo "count+compare launcher: $JG -> logs/AGE_SY/compare_v3_v6.out (at end)"
echo
squeue -u "$USER"
