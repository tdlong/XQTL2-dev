#!/bin/bash
# submit.build_v5.sh — refresh the pipeline to the annotated-catalog caller, then
# build the AGE_SY v5 founder catalog. One command:
#   bash scripts_oneoffs/AGE_SY/catalog_v5/submit.build_v5.sh
#
# Build only (~1.5 h founder recall, parallel per-chromosome). build_catalog.sh
# submits its own array + gather + filter chain and prints the job IDs. New defaults
# come from the pipeline: --min-dp 10 (#15), --snpgap 25 (#21), --maxaf 0.03,
# --exempt-founders B5:chr2L. Writes process/AGE_SY_v5/Catalog/ (catalog.annot.tsv.gz
# is the re-cuttable source; catalog.stats.txt is the per-rule tally to eyeball).
set -euo pipefail

# get the cluster's pipeline onto the annotated-catalog caller (5843d33, #21)
git -C pipeline pull --ff-only origin main
echo "pipeline now at: $(git -C pipeline log --oneline -1)"

bash pipeline/scripts/build_catalog.sh \
  --founders pipeline/helpfiles/B_founders.bams.txt \
  --out      process/AGE_SY_v5/Catalog
