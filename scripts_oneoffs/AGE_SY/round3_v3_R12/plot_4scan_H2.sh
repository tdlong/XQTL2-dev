#!/bin/bash
#SBATCH --job-name=plot4scanH2
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=30:00
#SBATCH -o figures/plot_4scan_H2.%j.out
###############################################################################
# plot_4scan_H2.sh — the 4-scan overlay figure with Cutler H² on the y-axis.
# Companion to plot_4scan.sh, which plots Wald_log10p from the same scan files.
# Reads SCAN_DIR from the environment (default process/AGE_SY); output is
# figures/<basename(SCAN_DIR)>_4scan_CutlH2.png.
#
#   sbatch --export=ALL,SCAN_DIR=process/AGE_SY \
#          scripts_oneoffs/AGE_SY/round3_v3_R12/plot_4scan_H2.sh
###############################################################################
module load R/4.2.2
mkdir -p figures
Rscript scripts_oneoffs/AGE_SY/common/plotting/aging_AGE_SY_H2.R
