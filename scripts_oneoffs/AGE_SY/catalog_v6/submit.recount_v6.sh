#!/bin/bash
# submit.recount_v6.sh — after XQTL2 #27 (catalog_count gets -I). The v6 catalog is
# already clean; this just re-counts the 44 samples + standard merge + compare vs v3.
# No re-filter, no recall.  One command:
#   bash scripts_oneoffs/AGE_SY/catalog_v6/submit.recount_v6.sh
set -euo pipefail
CATDIR=process/AGE_SY_v6/Catalog
DIR=process/AGE_SY_v6
BAMLIST=helpfiles/AGE_SY/bam_list.v4.txt
V3DIR=process/AGE_SY_v3
SELF=scripts_oneoffs/AGE_SY/catalog_v6

git -C pipeline pull --ff-only origin main
echo "pipeline: $(git -C pipeline log --oneline -1)"
[[ -s "$CATDIR/catalog.tsv.gz" ]] || { echo "no catalog at $CATDIR/catalog.tsv.gz" >&2; exit 1; }
# sanity: catalog must be duplicate-free before we trust the merge
dups=$(zcat "$CATDIR/catalog.tsv.gz" | cut -f1,2 | sort | uniq -d | wc -l)
echo "catalog duplicate positions: $dups (must be 0)"; [[ "$dups" -eq 0 ]] || { echo "catalog has dups — stop" >&2; exit 1; }

bash "$SELF/launch_count_compare.sh" "$DIR" "$CATDIR" "$BAMLIST" "$V3DIR"
squeue -u "$USER"
