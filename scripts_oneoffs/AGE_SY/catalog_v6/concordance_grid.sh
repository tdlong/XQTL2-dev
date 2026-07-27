#!/bin/bash
#SBATCH --job-name=concordance_grid
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:40:00
#SBATCH -o logs/AGE_SY/concordance_grid.out
# STEP 4 wrapper: run the 12-combo concordance (v6 vs v3). 4 cores = 24G for the merge.
#   sbatch concordance_grid.sh <v6 Calls dir> <v3 dir> <snp_pass.tsv.gz>
set -uo pipefail
module load R/4.2.2 2>/dev/null || true
Rscript scripts_oneoffs/AGE_SY/catalog_v6/concordance_grid.R "$1" "$2" "$3"
echo "done."
