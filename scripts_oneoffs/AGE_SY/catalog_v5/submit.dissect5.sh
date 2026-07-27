#!/bin/bash
#SBATCH --job-name=dissect5
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:30:00
#SBATCH -o logs/AGE_SY/dissect5.out
# Dissect 5 CLEAN (founder red-flag-free) v4-drops-alt SNPs for Con_R5_F, chrX.
# Each: REF/catalog/observed alleles, flag ladder + real callers, per-read
# base/baseQ/MAPQ split, tview. One log -> logs/AGE_SY/dissect5.out
#   bash scripts_oneoffs/AGE_SY/catalog_v5/submit.dissect5.sh   (self-sbatch below)
set -uo pipefail
D=scripts_oneoffs/AGE_SY/catalog_v5/dissect_snp.sh
for POS in 6471332 12335553 10242268 10928272 5332918; do
  echo; echo "==================================================================="
  CHR=chrX POS=$POS BAM=data/bam/AGE_SY/Con_R5_F.bam bash "$D"
done
