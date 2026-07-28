#!/bin/bash
#SBATCH --job-name=subset_compare
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:40:00
#SBATCH -o logs/AGE_SY/compare_v3_v7.out
# No recall, no recount: subset the existing loose v6 RefAlt to the final filters
# (maxaf 3%, snpgap 20), then run the pipeline's compare_refalt_calls.R vs v3.
#   sbatch scripts_oneoffs/AGE_SY/catalog_v7/submit.subset_compare.sh
set -uo pipefail
module load R/4.2.2 2>/dev/null || true
V6=process/AGE_SY_v6/Calls
PASS=process/AGE_SY_v6/Catalog/snp_pass.tsv.gz
V7=process/AGE_SY_v7/Calls
V3=process/AGE_SY_v3
for f in "$V6/RefAlt.chrX.txt" "$PASS"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

echo "=== subset loose v6 RefAlt -> final (maxaf 3%, snpgap 20) ==="
Rscript scripts_oneoffs/AGE_SY/catalog_v7/subset_final.R "$V6" "$PASS" "$V7" 20
echo "=== compare v3 vs v7 (final) ==="
Rscript pipeline/scripts/compare_refalt_calls.R "$V3" "$V7"
echo "done."
