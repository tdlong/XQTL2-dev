#!/bin/bash
#SBATCH --job-name=mal_part2
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0:20:00
#SBATCH -o logs/malathion/part2_%j.out
#SBATCH -e logs/malathion/part2_%j.out

# rerun_2022_part2.sh — stage 2 of rerun_2022_from_reads.sh: from BAMs and a
# built catalog to the finished diagnostics. It only submits jobs (bash + sbatch),
# so it is light; it runs as a job itself so it can wait on the catalog and the
# alignments with a SLURM dependency. Can also be run by hand on the login node
# once process/malathion_2022/Catalog/catalog.tsv.gz and the eight BAMs exist:
#   bash scripts_oneoffs/malathion/rerun_2022_part2.sh
set -euo pipefail
trap 'echo "ERROR: $0 aborted at line $LINENO (exit $?)" >&2' ERR
[ -d scripts_oneoffs/malathion ] || { echo "ERROR: run from the XQTL2-dev root" >&2; exit 1; }

PREFIX=malathion_2022
DIR=process/malathion_2022
CATALOG=$DIR/Catalog
BAMLIST=helpfiles/malathion/bam_list_2022.txt
DESIGN=helpfiles/malathion/design_2022.txt
PARFILE=helpfiles/malathion/hap_params_2022.R
SMOOTH=250
ACCOUNT=tdlong_lab
PARTITION=standard

# ── guards: stage 1 must be complete ─────────────────────────────────────────
[ -f "$CATALOG/catalog.tsv.gz" ] || { echo "ERROR: no $CATALOG/catalog.tsv.gz yet" >&2; exit 1; }
while read -r b; do
  [ -z "$b" ] || [ -e "$b" ] || { echo "ERROR: BAM missing: $b" >&2; exit 1; }
done < "$BAMLIST"
if [ -e "$DIR/Scans/${PREFIX}_all" ]; then
  echo "ERROR: $DIR/Scans/${PREFIX}_all exists; remove it or change PREFIX rather than overwrite" >&2; exit 1
fi
mkdir -p logs/malathion
echo "$(date)  stage 2: catalog and 8 BAMs present"

# ── 4. count against the catalog ─────────────────────────────────────────────
JID_REFALT=$(bash pipeline/scripts/call_samples.sh --catalog "$CATALOG" --bamlist "$BAMLIST" \
    --dir "$DIR" -A $ACCOUNT -p $PARTITION | tail -1)
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
echo "$(date)  stage 2 submitted. The finish job pushes when done."
