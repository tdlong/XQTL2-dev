#!/bin/bash
#SBATCH --job-name=plot4scan
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=30:00
#SBATCH -o figures/plot_4scan.%j.out
###############################################################################
# plot_4scan.sh — make the 4-scan overlay figure as a SLURM job.
# Reads SCAN_DIR from the environment (default process/AGE_SY_v3); output is
# figures/<basename(SCAN_DIR)>_4scan.png.
#
#   sbatch --export=ALL,SCAN_DIR=process/AGE_SY_v3_size75k \
#          scripts_oneoffs/AGE_SY/round3_v3_R12/plot_4scan.sh
###############################################################################
module load R/4.2.2
mkdir -p figures
Rscript scripts_oneoffs/AGE_SY/common/plotting/aging_AGE_SY_v3.R
