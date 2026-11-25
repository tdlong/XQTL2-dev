#!/bin/bash
# run_2022.sh — the standard pipeline on the eight 2022 malathion pools.
#
# Three steps, in order, from the repo root on the login node:
#   1. sbatch scripts_oneoffs/malathion/fetch_sra.sh      (16 files, md5-checked;
#        rerun it until  ls data/raw/malathion_2022/*.md5ok | wc -l  says 16)
#   2. wait for the A-founder catalog: process/malathion_2022/Catalog/catalog.tsv.gz
#        (build once with:  bash pipeline/scripts/build_catalog.sh
#            --founders pipeline/helpfiles/A_founders.bams.txt --out process/malathion_2022/Catalog)
#   3. bash scripts_oneoffs/malathion/run_2022.sh         (this: fq2bam -> calls -> haps -> scan)
#
# This is pipeline/scripts/run_full_pipeline.sh with the malathion settings
# filled in, nothing more. It checks 1 and 2 and refuses otherwise.
# helpfiles/malathion_2022/bam_list.txt is pre-made (8 samples + A founders)
# because run_full_pipeline.sh would otherwise draft it from BAMs that do not
# exist yet at submit time.
set -euo pipefail
shopt -s nullglob
[ -d scripts_oneoffs/malathion ] || { echo "ERROR: run from the XQTL2-dev root" >&2; exit 1; }

ok=(data/raw/malathion_2022/*.md5ok)
[ ${#ok[@]} -eq 16 ] || { echo "ERROR: ${#ok[@]} / 16 reads verified; run step 1 (fetch_sra.sh) until 16" >&2; exit 1; }
[ -f process/malathion_2022/Catalog/catalog.tsv.gz ] || { echo "ERROR: no catalog yet (step 2)" >&2; exit 1; }
[ -f helpfiles/malathion_2022/bam_list.txt ] || { echo "ERROR: missing helpfiles/malathion_2022/bam_list.txt" >&2; exit 1; }

bash pipeline/scripts/run_full_pipeline.sh \
    --project   malathion_2022 \
    --barcodes  helpfiles/malathion/malathion_2022.barcodes.txt \
    --rawdir    data/raw/malathion_2022 \
    --bamdir    data/bam/malathion_2022 \
    --catalog   process/malathion_2022/Catalog \
    --parfile   helpfiles/malathion/hap_params_2022.R \
    --design    helpfiles/malathion/design_2022.txt \
    --scan      malathion_2022_all \
    --founders  A
