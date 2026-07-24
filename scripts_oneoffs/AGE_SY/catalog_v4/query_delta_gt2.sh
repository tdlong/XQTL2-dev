#!/bin/bash
#SBATCH --job-name=query_dgt2
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:20:00
#SBATCH -o logs/AGE_SY/query_delta_gt2_%j.out
# Cheap QUERY of the already-built delta_gt2 tail (no re-merge). 4 cores = 24G for
# loading the ~30M-row tail table.
set -uo pipefail
module load R/4.2.2 2>/dev/null || true
Rscript scripts_oneoffs/AGE_SY/catalog_v4/query_delta_gt2.R process/AGE_SY_v4/Calls
echo "done."
