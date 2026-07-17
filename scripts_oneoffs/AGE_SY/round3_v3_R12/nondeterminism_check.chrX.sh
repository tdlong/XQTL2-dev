#!/bin/bash
#SBATCH --job-name=nondet_check
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00
#SBATCH -o process/prove_flips/nondet_task_%a.out
###############################################################################
# nondeterminism_check.chrX.sh
#
# The tiled vs non-tiled callers disagree at MANY SNPs, mid-tile, nowhere near
# boundaries. That rules out a padding/edge effect. Two candidates remain:
#   (a) the caller is NONDETERMINISTIC — two runs of the identical command give
#       different SNP sets scattered everywhere;
#   (b) a global (not edge) dependence on window extent.
# They split on one test: run the EXACT SAME command TWICE and diff the FULL
# filtered SNP tables (not a handful of sites — the whole table, at scale).
#
#   many scattered differences between two identical runs -> NONDETERMINISTIC.
#       (=> the tiled/non-tiled disagreement was never about tiling; it's noise)
#   two identical runs are identical                       -> deterministic;
#       the window-extent effect is real and we characterize it next.
#
# 4 arms as a parallel array: d=1000 x2 (A,B), d=50000 x2 (A,B), same -r window.
# Production uses -d 1000, so that pair is the one that matters; d=50000 is a
# cross-check. A dependent merge diffs A vs B within each depth.
#
# Submit once from repo root on the cluster (just bash it; it self-submits):
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/nondeterminism_check.chrX.sh
# Read the answer in process/prove_flips/nondet_merge.out
###############################################################################
set -uo pipefail

REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
REGION=chrX:15900000-16400000        # 500 kb contiguous -r window (~7k records)
OUT=process/prove_flips
# task index -> "depth rep"
C1="1000 A"; C2="1000 B"; C3="50000 A"; C4="50000 B"

# ======================= orchestrator (login node) ===========================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" && "${1:-}" != "--merge" ]]; then
  SELF=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
  mkdir -p "$OUT"
  JID=$(sbatch --parsable --array=1-4 "$SELF")
  echo "4 arms (d1000 A/B, d50000 A/B) submitted as a parallel array: $JID"
  MID=$(sbatch --parsable --dependency=afterok:${JID} \
        -A tdlong_lab -p standard --time=20:00 -o "$OUT/nondet_merge.out" \
        --wrap="bash '$SELF' --merge")
  echo "merge submitted: $MID"
  echo "watch:  squeue -u \$USER"
  echo "answer: $OUT/nondet_merge.out"
  exit 0
fi

# ============================ merge (verdict) ================================
if [[ "${1:-}" == "--merge" ]]; then
  echo "==================================================================="
  echo "NONDETERMINISM TEST — identical command run twice, full SNP tables"
  echo "REGION $REGION   filter -m2 -M2 -v snps -i 'QUAL>59'"
  echo "==================================================================="
  posfile() { cut -f1 "$1" | sort; }
  report() { # $1=A file  $2=B file  $3=label
    local A="$1" B="$2" lab="$3"
    local na nb aonly bonly both
    na=$(wc -l < "$A" | tr -d ' '); nb=$(wc -l < "$B" | tr -d ' ')
    aonly=$(comm -23 <(posfile "$A") <(posfile "$B") | wc -l | tr -d ' ')
    bonly=$(comm -13 <(posfile "$A") <(posfile "$B") | wc -l | tr -d ' ')
    both=$(comm -12 <(posfile "$A") <(posfile "$B") | wc -l | tr -d ' ')
    echo
    echo "----- $lab : two IDENTICAL runs -----"
    printf "  run A kept: %s SNPs\n  run B kept: %s SNPs\n" "$na" "$nb"
    printf "  in BOTH: %s   |   A-only: %s   B-only: %s   => DIFFER at %s positions\n" \
           "$both" "$aonly" "$bonly" "$((aonly+bonly))"
    if [[ "$((aonly+bonly))" -eq 0 ]]; then
      echo "  VERDICT: byte-identical -> DETERMINISTIC at this depth."
    else
      echo "  VERDICT: identical command, different SNP set -> NONDETERMINISTIC."
      echo "  example differing positions (QUAL in run A / run B; '-' = not called/failed):"
      comm -3 <(posfile "$A") <(posfile "$B") | tr -d '\t' | head -12 | while read -r p; do
        qa=$(awk -v p="$p" '$1==p{print $2}' "$A"); qb=$(awk -v p="$p" '$1==p{print $2}' "$B")
        printf "    chrX:%s   A=%-10s B=%-10s\n" "$p" "${qa:--}" "${qb:--}"
      done
    fi
  }
  report "$OUT/cond.1.txt" "$OUT/cond.2.txt" "d=1000 (production depth)"
  report "$OUT/cond.3.txt" "$OUT/cond.4.txt" "d=50000"
  echo
  echo "done."
  exit 0
fi

# ============================ array task (one arm) ==========================
module load bcftools/1.21 2>/dev/null || true
eval "cond=\$C$SLURM_ARRAY_TASK_ID"
read -r depth rep <<< "$cond"
vcf="$OUT/nondet.d${depth}.${rep}.vcf.gz"
echo "[$(date +%H:%M:%S)] arm d$depth $rep  $REGION" >&2
bcftools mpileup -I -d "$depth" -r "$REGION" -a FORMAT/AD,FORMAT/DP -f "$REF" -b "$BAMS" \
  | bcftools call -mv -Oz -o "$vcf"
bcftools index -f "$vcf"
n=$(bcftools view -H "$vcf" 2>/dev/null | wc -l | tr -d ' ')
[[ "$n" -eq 0 ]] && { echo "ERROR: 0 records (tool failure) for d$depth $rep" >&2; exit 1; }
# passing SNPs: POS<TAB>QUAL, sorted by POS
bcftools view -m2 -M2 -v snps -i 'QUAL>59' "$vcf" 2>/dev/null \
  | bcftools query -f '%POS\t%QUAL\n' | sort -k1,1n > "$OUT/cond.$SLURM_ARRAY_TASK_ID.txt"
echo "[$(date +%H:%M:%S)] arm d$depth $rep done ($(wc -l < "$OUT/cond.$SLURM_ARRAY_TASK_ID.txt") SNPs)" >&2
