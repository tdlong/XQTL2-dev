#!/bin/bash
#SBATCH --job-name=merge_concord
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=6G
#SBATCH --time=01:00:00
#SBATCH -o logs/AGE_SY/concordance_grid.out
# Finish the pipeline after the count-merge crash (XQTL2 #24): dedup-merge the existing
# v6 counts into RefAlt.<chr>.txt (dropping multiallelic positions), then run the
# 12-combo concordance vs v3.  ->  logs/AGE_SY/concordance_grid.out
#   sbatch scripts_oneoffs/AGE_SY/catalog_v6/submit.merge_concordance.sh
set -uo pipefail
module load R/4.2.2 2>/dev/null || true
CALLS=process/AGE_SY_v6/Calls
PASS=process/AGE_SY_v6/Catalog/snp_pass.tsv.gz
V3=process/AGE_SY_v3
echo "=== dedup-merge v6 counts -> RefAlt ==="
Rscript scripts_oneoffs/AGE_SY/catalog_v6/merge_dedup.R "$CALLS"
echo "=== concordance v6 vs v3 (12 combos) ==="
Rscript scripts_oneoffs/AGE_SY/catalog_v6/concordance_grid.R "$CALLS" "$V3" "$PASS"
echo "done."
