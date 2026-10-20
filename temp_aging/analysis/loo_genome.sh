#!/bin/bash
#SBATCH --job-name=loo_genome
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=4:00:00
#SBATCH -o logs/loo_genome.out
#SBATCH -e logs/loo_genome.err

# Submit from the XQTL2-dev root:  sbatch temp_aging/analysis/loo_genome.sh
# 24 GB total. Each meansBySample is ~257 MB on disk and is read one at a time,
# reduced, then dropped -- the four are never in memory together.

module load R/4.2.2
Rscript temp_aging/analysis/loo_genome.R
