#!/bin/bash
#SBATCH --job-name=agree_delta
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:40:00
#SBATCH -o logs/AGE_SY/agreement_delta_%j.out
# 4 cores = 24G (memory, not compute): merges the two RefAlt tables per chromosome.
#
# Deterministic per-(SNP,sample) count agreement, v3 vs v4/Calls, on the shared
# SNPs and shared samples. Same BAM -> counts should be identical; this reports
# what fraction are exactly equal and, where not, the delta distribution vs depth.
#
#   sbatch scripts_oneoffs/AGE_SY/catalog_v4/agreement_delta.sh
set -uo pipefail
module load R/4.2.2 2>/dev/null || true
Rscript scripts_oneoffs/AGE_SY/catalog_v4/agreement_delta.R \
    process/AGE_SY_v3 process/AGE_SY_v4/Calls
echo "done."
