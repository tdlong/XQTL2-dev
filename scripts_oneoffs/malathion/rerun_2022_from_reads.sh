#!/bin/bash
###############################################################################
# rerun_2022_from_reads.sh — the full 2022 malathion experiment, all 8 pools,
# from the public reads to the diagnostics, through the current pipeline.
#
# Rep 2 was never on HPC3 (logs/malathion/inventory.txt). All eight pools are on
# the SRA under BioProject PRJNA857080, so this fetches them and runs every step
# of pipeline/README.md with SLURM dependencies, the way run_full_pipeline.sh
# chains them, plus the per-pair scans and the paper's ANOVA from rerun_2022.sh:
#
#   stage 1 (this script, login node; everything below is sbatch)
#   fetch_sra.sh          8 runs x 2 mates from ENA, md5  -> data/raw/malathion_2022
#   build_catalog.sh      A-founder catalog (none on HPC3; AGE_SY's is B founders)
#   fq2bam.sh             align, one job per sample       -> data/bam/malathion_2022
#   stage 2 (rerun_2022_part2.sh, an sbatch job that waits on stage 1)
#   call_samples.sh       count against the catalog       -> process/malathion_2022/Calls
#   refalt_qc / bam_qc    advisory, not in the chain      -> logs/malathion/qc_*.out
#   run_haps.sh           haplotypes, 8 pools             -> process/malathion_2022/Haps
#   run_scan.sh x5        all pools + one scan per pair   -> process/malathion_2022/Scans
#   finish_2022.sh        paper ANOVA, Fig.4 tables       -> logs/malathion/, git push
#
# Two stages because call_samples.sh checks that catalog.tsv.gz exists when it
# is SUBMITTED, so it cannot be chained behind a catalog that is still building.
# Stage 2 is therefore submitted as a job that starts after the catalog and the
# alignments have finished, and submits the rest from there.
#
# Everything is relative to the repo root; nothing touches an existing folder.
# Verified reads, an existing catalog and existing BAMs are reused; anything
# downstream is recomputed. Do not rerun while its jobs are still in the queue.
#
# Login node: bash and sbatch only. Run once from the repo root:
#   bash scripts_oneoffs/malathion/rerun_2022_from_reads.sh
#
# Wall time: download and catalog in parallel (hours), alignment up to a day per
# sample in parallel, then calls, haplotypes (up to a day on highmem) and
# scans. Two to three days end to end.
###############################################################################
set -euo pipefail
trap 'echo "ERROR: $0 aborted at line $LINENO (exit $?)" >&2' ERR
shopt -s nullglob
[ -d scripts_oneoffs/malathion ] || { echo "ERROR: run from the XQTL2-dev root" >&2; exit 1; }

PROJECT=malathion_2022
PREFIX=malathion_2022
DIR=process/$PROJECT
RAW=data/raw/$PROJECT
BAM=data/bam/$PROJECT
CATALOG=$DIR/Catalog
FOUNDERS=pipeline/helpfiles/A_founders.bams.txt
SRA_TSV=helpfiles/malathion/sra_runs.tsv
BARCODES=helpfiles/malathion/malathion_2022.barcodes.txt
BAMLIST=helpfiles/malathion/bam_list_2022.txt
DESIGN=helpfiles/malathion/design_2022.txt
PARFILE=helpfiles/malathion/hap_params_2022.R
SMOOTH=250
ACCOUNT=tdlong_lab
PARTITION=standard

# ── guards ───────────────────────────────────────────────────────────────────
for f in "$SRA_TSV" "$BARCODES" "$BAMLIST" "$DESIGN" "$PARFILE" "$FOUNDERS"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done
for r in 1 2 3 4; do
  [ -f "helpfiles/malathion/design_2022_R$r.txt" ] || { echo "ERROR: missing design_2022_R$r.txt" >&2; exit 1; }
done
[ -e pipeline/ref/dm6.fa ] || { echo "ERROR: pipeline/ref/dm6.fa not found (fq2bam.sh needs it)" >&2; exit 1; }
while read -r b; do
  [ -z "$b" ] || [ -e "$b" ] || { echo "ERROR: founder BAM missing: $b" >&2; exit 1; }
done < "$FOUNDERS"
NRUN=$(($(wc -l < "$SRA_TSV") - 1))
NBC=$(wc -l < "$BARCODES")
[ "$NRUN" -eq 8 ] && [ "$NBC" -eq 8 ] || { echo "ERROR: expected 8 runs and 8 barcode rows, got $NRUN and $NBC" >&2; exit 1; }
if [ -e "$DIR/Scans/${PREFIX}_all" ]; then
  echo "ERROR: $DIR/Scans/${PREFIX}_all exists; remove it or change PREFIX rather than overwrite" >&2; exit 1
fi
mkdir -p "$RAW" "$BAM" "$DIR" logs/malathion

# ── 1. reads: a batch array, one run per task, both mates, md5-checked ────────
# A task exits non-zero unless both of its files match ENA's md5, and every job
# below waits on the whole array with afterok, so nothing can start on a partial
# or corrupt file. A verified file carries a .md5ok marker; with all 16 present
# the fetch is skipped. Never run two fetches into the same folder at once.
ok_files=("$RAW"/*.md5ok); n_ok=${#ok_files[@]}
if [ "$n_ok" -eq 16 ]; then
  JID_FETCH=""
  echo "reads:      16 / 16 verified in $RAW, fetch skipped"
else
  if squeue -u "$USER" -h -n sra_fetch | grep -q .; then
    echo "ERROR: an sra_fetch job is already in the queue; wait for it or scancel it first" >&2; exit 1
  fi
  JID_FETCH=$(sbatch --parsable -A $ACCOUNT scripts_oneoffs/malathion/fetch_sra.sh | cut -d_ -f1)
  echo "fetch:      $JID_FETCH  (8 runs x 2 mates from ENA -> $RAW, md5-checked; $n_ok / 16 already verified)"
fi

# ── 2. catalog, in parallel with the download ────────────────────────────────
if [ -f "$CATALOG/catalog.tsv.gz" ]; then
  JID_CAT=""
  echo "catalog:    exists, reused ($CATALOG)"
elif [ -d "$CATALOG/work" ]; then
  echo "ERROR: $CATALOG is still being built (work/ present, no catalog.tsv.gz)." >&2
  echo "       Wait for the catalog jobs to finish (squeue -u \$USER), then rerun." >&2
  echo "       If they failed, remove $CATALOG and rerun to rebuild." >&2
  exit 1
else
  JID_CAT=$(bash pipeline/scripts/build_catalog.sh --founders "$FOUNDERS" --out "$CATALOG" -A $ACCOUNT -p $PARTITION | tail -1)
  echo "catalog:    $JID_CAT  (A founders -> $CATALOG)"
fi

# ── 3. align ─────────────────────────────────────────────────────────────────
nbam=0
while read -r b; do [ -e "$b" ] && nbam=$((nbam+1)); done < <(head -8 "$BAMLIST")
if [ "$nbam" -eq 8 ]; then
  JID_BAM=""
  echo "fq2bam:     all 8 BAMs exist in $BAM, reused"
else
  JID_BAM=$(sbatch --parsable ${JID_FETCH:+--dependency=afterok:$JID_FETCH} -A $ACCOUNT \
      --array=1-$NBC pipeline/scripts/fq2bam.sh "$BARCODES" "$RAW" "$BAM" | cut -d_ -f1)
  echo "fq2bam:     $JID_BAM  (8 samples -> $BAM)"
fi

# ── stage 2, once the catalog and the BAMs are there ────────────────────────
DEP=$(printf '%s:%s' "$JID_CAT" "$JID_BAM" | sed 's/^://; s/:$//')
JID_P2=$(sbatch --parsable ${DEP:+--dependency=afterok:$DEP} -A $ACCOUNT \
    scripts_oneoffs/malathion/rerun_2022_part2.sh)
echo "stage 2:    $JID_P2  (waits on ${DEP:-nothing}; submits calls, haplotypes, scans, finish)"
echo
echo "Stage 1 submitted. Stage 2 submits the rest when the catalog and BAMs exist;"
echo "its log is logs/malathion/part2_$JID_P2.out. The finish job pushes when all is done."

