#!/bin/bash
#SBATCH --job-name=v3_eff_compare
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:30:00
#SBATCH -o logs/AGE_SY/v3_effective_compare_%j.out
###############################################################################
# v3_effective_compare.sh — three-way SNP-set comparison, needs no v4 calls:
#   v3_raw (QUAL>59)  vs  v3_effective (v3 after good_SNPs)  vs  catalog (v4).
# The "effective" v3 is what the haplotyper actually sees; the catalog should
# reproduce it. Runs on existing files (v3 RefAlt + the built catalog).
#
#   sbatch scripts_oneoffs/AGE_SY/catalog_v4/v3_effective_compare.sh
###############################################################################
set -uo pipefail
module load R/4.2.2 2>/dev/null || true

Rscript scripts_oneoffs/AGE_SY/catalog_v4/v3_effective_compare.R \
    process/AGE_SY_v3 \
    process/AGE_SY_v4/Catalog/catalog.tsv.gz \
    helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R
echo "done."
