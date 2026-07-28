#!/bin/bash
#SBATCH --job-name=compare_summary
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:40:00
#SBATCH -o logs/AGE_SY/compare_summary.out
# Better v3-vs-v6 summary: overlap counts + |freq diff| joint distribution vs coverage
# (+ a hexbin plot). Reads the clean RefAlt already on disk.
#   sbatch scripts_oneoffs/AGE_SY/catalog_v6/submit.compare_summary.sh
set -uo pipefail
module load R/4.4.2 2>/dev/null || true
Rscript scripts_oneoffs/AGE_SY/catalog_v6/compare_summary.R process/AGE_SY_v3 process/AGE_SY_v6/Calls
echo "done."
