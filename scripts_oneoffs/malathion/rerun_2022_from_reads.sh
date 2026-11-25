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
#   fetch_sra.sh          8 runs from the ENA mirror     -> data/raw/malathion_2022
#   build_catalog.sh      A-founder catalog (none on HPC3; AGE_SY's is B founders)
#   fq2bam.sh             align, one job per sample       -> data/bam/malathion_2022
#   call_samples.sh       count against the catalog       -> process/malathion_2022/Calls
#   refalt_qc / bam_qc    advisory, not in the chain      -> logs/malathion/qc_*.out
#   run_haps.sh           haplotypes, 8 pools             -> process/malathion_2022/Haps
#   run_scan.sh x5        all pools + one scan per pair   -> process/malathion_2022/Scans
#   finish_2022.sh        paper ANOVA, Fig.4 tables       -> logs/malathion/, git push
#
# Everything is relative to the repo root; nothing touches an existing folder.
# Idempotent where it is cheap: fetched files with a good md5, an existing
# catalog, and existing BAMs are reused; anything downstream is recomputed.
#
# Login node: bash and sbatch only. Run once from the repo root:
#   bash scripts_oneoffs/malathion/rerun_2022_from_reads.sh
#
# Wall time: the download and the catalog run in parallel (hours), alignment is
# up to a day per sample in parallel, then calls, haplotypes (up to a day on
# highmem) and scans. Expect two to three days end to end.
###############################################################################
set -euo pipefail
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

# ── 1. reads ─────────────────────────────────────────────────────────────────
JID_FETCH=$(sbatch --parsable -A $ACCOUNT scripts_oneoffs/malathion/fetch_sra.sh)
echo "fetch:      $JID_FETCH  (8 runs, ENA -> $RAW)"

# ── 2. catalog, in parallel with the download ────────────────────────────────
if [ -f "$CATALOG/catalog.tsv.gz" ]; then
  JID_CAT=""
  echo "catalog:    exists, reused ($CATALOG)"
else
  JID_CAT=$(bash pipeline/scripts/build_catalog.sh --founders "$FOUNDERS" --out "$CATALOG" -A $ACCOUNT -p $PARTITION | tail -1)
  echo "catalog:    $JID_CAT  (A founders -> $CATALOG)"
fi

# ── 3. align ─────────────────────────────────────────────────────────────────
JID_BAM=$(sbatch --parsable --dependency=afterok:$JID_FETCH -A $ACCOUNT \
    --array=1-$NBC pipeline/scripts/fq2bam.sh "$BARCODES" "$RAW" "$BAM" | cut -d_ -f1)
echo "fq2bam:     $JID_BAM  (8 samples -> $BAM)"

# ── 4. count against the catalog ─────────────────────────────────────────────
DEP_CALL="$JID_BAM${JID_CAT:+:$JID_CAT}"
JID_REFALT=$(bash pipeline/scripts/call_samples.sh --catalog "$CATALOG" --bamlist "$BAMLIST" \
    --dir "$DIR" --after "$DEP_CALL" -A $ACCOUNT -p $PARTITION | tail -1)
echo "RefAlt:     $JID_REFALT  (-> $DIR/Calls)"

# advisory QC, as run_full_pipeline.sh does; never gates the chain
QC="module load R/4.2.2 && Rscript pipeline/scripts/refalt_qc.R --dir $DIR --parfile $PARFILE"
QC="$QC && module load samtools/1.10 && bash pipeline/scripts/bam_qc.sh --bamlist $BAMLIST --out $DIR/Calls/bam_qc.txt"
JID_QC=$(sbatch --parsable --dependency=afterok:$JID_REFALT -A $ACCOUNT -p $PARTITION \
    --cpus-per-task=1 --mem-per-cpu=8G --time=2:00:00 -o logs/malathion/qc_%j.out --wrap="$QC")
echo "sample QC:  $JID_QC  (-> $DIR/Calls/refalt_qc.txt, logs/malathion/qc_$JID_QC.out)"

# ── 5. haplotypes ────────────────────────────────────────────────────────────
JID_HAPS=$(bash pipeline/scripts/run_haps.sh --parfile "$PARFILE" --dir "$DIR" --after "$JID_REFALT" -A $ACCOUNT)
echo "haplotypes: $JID_HAPS  (REFALT2haps, 5 chromosomes, highmem)"

# ── 6. scans: all pools, then one per pair (same as rerun_2022.sh) ──────────
submit_scan () {  # submit_scan <name> <design> <sex>
  local out jid
  out=$(bash pipeline/scripts/run_scan.sh --design "$2" --dir "$DIR" --scan "$1" \
          --smooth "$SMOOTH" --sex "$3" --after "$JID_HAPS")
  jid=$(printf '%s\n' "$out" | awk '/^done:/{print $2}')
  [ -n "$jid" ] || { echo "ERROR: no concat job id for $1" >&2; printf '%s\n' "$out" >&2; exit 1; }
  JIDS="${JIDS:+$JIDS:}$jid"
  echo "scan:       $1  ($2, --sex $3, concat job $jid)"
}
JIDS=""; ALSO=""
submit_scan "${PREFIX}_all" "$DESIGN" mixed
for r in 1 2 3 4; do
  dr="helpfiles/malathion/design_2022_R$r.txt"
  sex=mixed
  if   grep -q '_F_' "$dr" && ! grep -q '_M_' "$dr"; then sex=F
  elif grep -q '_M_' "$dr" && ! grep -q '_F_' "$dr"; then sex=M; fi
  submit_scan "${PREFIX}_R$r" "$dr" "$sex"
  ALSO="${ALSO:+$ALSO,}${PREFIX}_R$r"
done

# ── 7. diagnostics + return through git ──────────────────────────────────────
JID_FIN=$(sbatch --parsable --dependency=afterok:"$JIDS" -A $ACCOUNT \
    scripts_oneoffs/malathion/finish_2022.sh "$DIR" "$PREFIX" "$ALSO")
echo "finish:     $JID_FIN  (paper ANOVA, peak tables, figures -> logs/malathion/, then git push)"
echo
echo "Everything is submitted. squeue -u \$USER shows the chain; the finish job pushes when done."
