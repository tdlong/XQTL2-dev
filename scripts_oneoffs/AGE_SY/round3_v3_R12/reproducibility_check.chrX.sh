#!/bin/bash
#SBATCH --job-name=repro_qual
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=01:00:00
#SBATCH -o process/prove_flips/repro_task_%a.out
###############################################################################
# reproducibility_check.chrX.sh  (fast: parallel array + narrow windows)
#
# prove_depth_flips found 5/7 flips are the QUAL>59 gate, and at 3 sites the
# per-sample DP:AD is IDENTICAL at d=1000 vs d=50000 yet QUAL crosses 59. This
# decides, with no speculation, which is true:
#   (1) NONDETERMINISTIC: QUAL differs between two runs at the SAME -d
#       => the flips are run-to-run noise, not a real depth effect.
#   (2) DETERMINISTIC: two same--d runs give identical QUAL
#       => real -d effect; the full INFO records (dumped below) show which
#          field carries the change.
#
# Speed: a call at a position depends only on reads OVERLAPPING it, so instead
# of the 500 kb window we call ±10 kb around each of the 7 sites (~60 kb total)
# and run the 4 arms (d1000 x2, d50000 x2) as a PARALLEL array. The merge prints
# each site's QUAL next to the value the 500 kb run gave (EXPECTED_*), so if the
# narrow window changed any call it is caught, not hidden.
#
# Submit once from repo root on the cluster (just bash it; it self-submits):
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/reproducibility_check.chrX.sh
# Read the answer in process/prove_flips/repro_merge.out
###############################################################################
set -uo pipefail

REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
OUT=process/prove_flips
POS="16025382 16030427 16030958 16031180 16031851 16045058 16363526"
# ±10 kb windows around the 7 sites, merged to non-overlapping regions:
REGIONS="chrX:16015382-16055058,chrX:16353526-16373526"
# task index -> "depth tag"
C1="1000 d1000a"; C2="1000 d1000b"; C3="50000 d50000a"; C4="50000 d50000b"
# QUAL the earlier ±10 kb run (depth_cap_check) produced — SAME window as this
# script, so if the caller is deterministic the a/b columns must equal these.
# (They do NOT match the 500 kb run: window changes QUAL — that's a separate,
# already-established finding.)
EXPECTED_d1000="16025382=36.3391 16030427=96.212 16030958=390.393 16031180=110.727 16031851=134.905 16045058=64.0468 16363526=24.2829"
EXPECTED_d50000="16025382=66.206 16030427=132.212 16030958=376.591 16031180=85.5538 16031851=111.443 16045058=67.0423 16363526=36.8386"

# ======================= orchestrator (login node) ===========================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" && "${1:-}" != "--merge" ]]; then
  SELF=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
  mkdir -p "$OUT"
  JID=$(sbatch --parsable --array=1-4 "$SELF")
  echo "4 arms submitted as a parallel array: $JID"
  MID=$(sbatch --parsable --dependency=afterok:${JID} \
        -A tdlong_lab -p standard --time=15:00 -o "$OUT/repro_merge.out" \
        --wrap="bash '$SELF' --merge")
  echo "merge submitted: $MID"
  echo "watch:  squeue -u \$USER"
  echo "answer: $OUT/repro_merge.out"
  exit 0
fi

# ============================ merge (verdict) ================================
if [[ "${1:-}" == "--merge" ]]; then
  module load bcftools/1.21 2>/dev/null || true
  echo "==================================================================="
  echo "QUAL at the 7 flip sites — 4 runs (narrow ±10 kb windows)"
  echo "  same-depth cols EQUAL  => deterministic (real -d effect)"
  echo "  same-depth cols DIFFER => nondeterministic (flips are RNG noise)"
  echo "  EXPECT cols = what the 500 kb run gave; must match (validates the shrink)"
  echo "==================================================================="
  printf "%-10s %-11s %-11s | %-11s %-11s | %-9s %-9s\n" \
         "POS" "d1000a" "d1000b" "d50000a" "d50000b" "EXP1000" "EXP50000"
  get() { awk -v p="$2" '$1==p{print ($3==""?"no-call":$3)}' "$OUT/repro.$1.persite.tsv" 2>/dev/null | head -1; }
  exp() { echo "$1" | tr ' ' '\n' | awk -F= -v p="$2" '$1==p{print $2}'; }
  for p in $POS; do
    printf "%-10s %-11s %-11s | %-11s %-11s | %-9s %-9s\n" "$p" \
      "$(get d1000a "$p")" "$(get d1000b "$p")" "$(get d50000a "$p")" "$(get d50000b "$p")" \
      "$(exp "$EXPECTED_d1000" "$p")" "$(exp "$EXPECTED_d50000" "$p")"
  done
  echo
  echo "Full site record (VCF cols 1-8, INFO included) d=1000a vs d=50000a —"
  echo "which INFO field carries the QUAL swing at the identical-count sites:"
  for p in $POS; do
    echo "--- chrX:$p ---"
    printf "  d1000 : %s\n" "$(awk -v p="$p" '$1==p{sub(/^[0-9]+\t/,"");print}' "$OUT/repro.d1000a.rec.tsv" 2>/dev/null)"
    printf "  d50000: %s\n" "$(awk -v p="$p" '$1==p{sub(/^[0-9]+\t/,"");print}' "$OUT/repro.d50000a.rec.tsv" 2>/dev/null)"
  done
  echo
  echo "done."
  exit 0
fi

# ============================ array task (one arm) ==========================
module load bcftools/1.21 2>/dev/null || true
eval "cond=\$C$SLURM_ARRAY_TASK_ID"
read -r depth tag <<< "$cond"
vcf="$OUT/repro.$tag.vcf.gz"
echo "[$(date +%H:%M:%S)] arm $tag  -d $depth  $REGIONS" >&2
bcftools mpileup -I -d "$depth" -r "$REGIONS" -a FORMAT/AD,FORMAT/DP -f "$REF" -b "$BAMS" \
  | bcftools call -mv -Oz -o "$vcf"
bcftools index -f "$vcf"

: > "$OUT/repro.$tag.persite.tsv"
: > "$OUT/repro.$tag.rec.tsv"
for p in $POS; do
  # per-site: POS  nALT  QUAL   (blank QUAL => not called)
  bcftools query -r "chrX:$p-$p" -f '%POS\t%ALT\t%QUAL\n' "$vcf" 2>/dev/null \
    | awk 'NR==1{na=($2=="."?0:split($2,a,",")); print $1"\t"na"\t"$3}' >> "$OUT/repro.$tag.persite.tsv"
  # full record (cols 1-8) prefixed by POS for lookup
  rec=$(bcftools view -H -r "chrX:$p-$p" "$vcf" 2>/dev/null | cut -f1-8 | head -1)
  printf "%s\t%s\n" "$p" "$rec" >> "$OUT/repro.$tag.rec.tsv"
done
echo "[$(date +%H:%M:%S)] arm $tag done" >&2
