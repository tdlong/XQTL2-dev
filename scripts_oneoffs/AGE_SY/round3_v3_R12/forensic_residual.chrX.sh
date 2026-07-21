#!/bin/bash
#SBATCH --job-name=forensic_resid
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00
#SBATCH -o process/prove_flips/forensic_%a.out
###############################################################################
# forensic_residual.chrX.sh
#
# THE unexplained thing: at 16045058 (uncapped, every sample far under -d 1000)
# the whole-chr caller gives QUAL 58.07 and the tiled caller keeps it (>59) —
# same reads, very different QUAL. "The window changes QUAL" is a restatement,
# not a cause. QUAL comes from per-sample genotype likelihoods, which depend on
# BASE QUALITIES, and BAQ (base alignment quality, ON by default) is the one
# thing that rewrites base qualities. Hypothesis: BAQ is computed window-
# dependently -> same reads, different base quals, different QUAL.
#
# Call 16045058 across window sizes, BAQ on AND off, dump the FULL record each
# time (QUAL + per-sample AD + PL). We then read, directly:
#   (1) is the per-sample AD byte-identical across windows?  (is the data same?)
#   (2) does QUAL actually move with window size?            (reproduce it)
#   (3) does -B (BAQ off) remove the window-dependence?      (is BAQ the cause?)
#
# 6 conditions {50k,250k,1M half-width} x {BAQ on, off} as a parallel array + a
# dependent merge. (5 Mb dropped: 0.1/0.5/2 Mb already span the far-from-boundary
# range and it doesn't add anything.) Whole-chr reference: calls.chrX.bcf (58.07).
#
# Submit from repo root (self-submits):
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/forensic_residual.chrX.sh
# Read: process/prove_flips/forensic_merge.out
###############################################################################
set -uo pipefail

REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
BCF=process/AGE_SY_v3/calls.chrX.bcf     # whole-chr reference record (kept)
OUT=process/prove_flips
POS=16045058
# condition index -> "halfwidth baqflag label"
C1="50000  '' BAQon_50k";    C2="50000  -B BAQoff_50k"
C3="250000 '' BAQon_250k";   C4="250000 -B BAQoff_250k"
C5="1000000 '' BAQon_1M";    C6="1000000 -B BAQoff_1M"

# ======================= orchestrator =======================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" && "${1:-}" != "--merge" ]]; then
  SELF=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
  mkdir -p "$OUT"
  JID=$(sbatch --parsable --array=1-6 "$SELF")
  echo "6 conditions submitted as a parallel array: $JID"
  MID=$(sbatch --parsable --dependency=afterok:${JID} \
        -A tdlong_lab -p standard --time=20:00 -o "$OUT/forensic_merge.out" \
        --wrap="bash '$SELF' --merge")
  echo "merge submitted: $MID   -> $OUT/forensic_merge.out"
  exit 0
fi

# ======================= merge =======================
if [[ "${1:-}" == "--merge" ]]; then
  module load bcftools/1.21 2>/dev/null || true
  echo "==================================================================="
  echo "FORENSIC on chrX:$POS   (residual: uncapped, QUAL flips with window)"
  echo "==================================================================="
  echo "whole-chr reference (calls.chrX.bcf):"
  bcftools query -r "chrX:$POS-$POS" -f '  QUAL=%QUAL  ALT=%ALT  INFO/DP=%INFO/DP\n' "$BCF" 2>/dev/null | head -1
  echo
  printf "%-14s %-8s %-10s %-6s  %s\n" "condition" "QUAL" "nALT" "ADeq?" "(ADeq = per-sample AD identical to BAQon_1M)"
  ref="$OUT/forensic.ad.5.txt"   # BAQon_1M (largest window) as the reference AD vector
  for i in $(seq 1 6); do
    eval "c=\$C$i"; read -r w b lab <<< "$c"
    q=$(cat "$OUT/forensic.q.$i.txt" 2>/dev/null); na=$(cat "$OUT/forensic.na.$i.txt" 2>/dev/null)
    adeq="?"; if [[ -f "$OUT/forensic.ad.$i.txt" && -f "$ref" ]]; then
      cmp -s "$OUT/forensic.ad.$i.txt" "$ref" && adeq="YES" || adeq="no"; fi
    printf "%-14s %-8s %-10s %-6s\n" "$lab" "${q:--}" "${na:--}" "$adeq"
  done
  echo
  echo "read this:"
  echo "  * ADeq=YES across all rows  -> the reads/counts are the SAME regardless of window."
  echo "  * QUAL still varies down the BAQon column -> window changes QUAL on identical data."
  echo "  * BAQoff QUAL constant across window sizes -> BAQ is the cause."
  echo
  echo "full records per condition (QUAL/INFO/first-3-sample FORMAT) in $OUT/forensic_<a>.out"
  exit 0
fi

# ======================= array task =======================
module load bcftools/1.21 2>/dev/null || true
eval "c=\$C$SLURM_ARRAY_TASK_ID"; read -r W BAQ LAB <<< "$c"
[[ "$BAQ" == "''" ]] && BAQ=""     # empty flag
rs=$((POS - W)); [[ $rs -lt 1 ]] && rs=1; re=$((POS + W))
echo "[$(date +%H:%M:%S)] cond $SLURM_ARRAY_TASK_ID $LAB  -r chrX:$rs-$re  BAQ='${BAQ:-on}'" >&2
vcf="$OUT/forensic.$SLURM_ARRAY_TASK_ID.vcf.gz"
bcftools mpileup -I -d 1000 $BAQ -r "chrX:$rs-$re" -a FORMAT/AD,FORMAT/DP -f "$REF" -b "$BAMS" \
  | bcftools call -mv -Oz -o "$vcf"
bcftools index -f "$vcf"
# per-sample AD vector (to test data-identity), QUAL, nALT
bcftools query -r "chrX:$POS-$POS" -f '[%AD ]\n'  "$vcf" 2>/dev/null | head -1 > "$OUT/forensic.ad.$SLURM_ARRAY_TASK_ID.txt"
bcftools query -r "chrX:$POS-$POS" -f '%QUAL\n'   "$vcf" 2>/dev/null | head -1 > "$OUT/forensic.q.$SLURM_ARRAY_TASK_ID.txt"
bcftools query -r "chrX:$POS-$POS" -f '%ALT\n'    "$vcf" 2>/dev/null | head -1 \
  | awk '{print ($1=="."?0:split($1,a,","))}' > "$OUT/forensic.na.$SLURM_ARRAY_TASK_ID.txt"
# full record (QUAL/INFO + first 3 samples) for manual diff
echo "== cond $SLURM_ARRAY_TASK_ID $LAB =="
bcftools view -H -r "chrX:$POS-$POS" "$vcf" 2>/dev/null | cut -f1-12
