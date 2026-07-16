#!/bin/bash
#SBATCH --job-name=repro_qual
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00
#SBATCH -o process/prove_flips/repro_%j.out
###############################################################################
# reproducibility_check.chrX.sh
#
# prove_depth_flips proved 5/7 flips are the QUAL>59 gate, and — the sharp part —
# at 16031851 / 16045058 / 16363526 the per-sample DP:AD is IDENTICAL at d=1000
# vs d=50000 yet QUAL swings past 59. This script decides, with no speculation,
# WHICH of two things is true:
#
#   (1) NONDETERMINISTIC: QUAL differs between two runs at the SAME -d.
#       => the "flips" are run-to-run noise, not a real depth effect.
#   (2) DETERMINISTIC: two runs at the same -d give identical QUAL.
#       => it's a real -d effect; the full INFO-record dump below shows which
#          field (DP, QS, I16, MQ...) carries the QUAL change.
#
# Runs the EXACT production caller 4x over the same -r window:
#   d=1000 twice (a,b), d=50000 twice (a,b).
#
# Submit from repo root on the cluster:
#   mkdir -p process/prove_flips
#   sbatch scripts_oneoffs/AGE_SY/round3_v3_R12/reproducibility_check.chrX.sh
###############################################################################
set -uo pipefail
# only bcftools/1.21 (htslib>=1.16); do NOT load samtools/1.10 alongside it.
module load bcftools/1.21 2>/dev/null || true

REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
REGION=chrX:15900000-16400000
OUT=process/prove_flips
mkdir -p "$OUT"
POS="16025382 16030427 16030958 16031180 16031851 16045058 16363526"

for f in "$REF" "$BAMS"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

run() { # $1=depth  $2=tag
  local vcf="$OUT/repro.$2.vcf.gz"
  echo "[$(date +%H:%M:%S)] run $2  (-d $1) ..." >&2
  bcftools mpileup -I -d "$1" -r "$REGION" -a FORMAT/AD,FORMAT/DP -f "$REF" -b "$BAMS" \
    | bcftools call -mv -Oz -o "$vcf"
  bcftools index -f "$vcf"
  local n; n=$(bcftools view -H "$vcf" 2>/dev/null | wc -l | tr -d ' ')
  echo "  -> $n records" >&2
  [[ "$n" -eq 0 ]] && { echo "ERROR: 0 records for $2 — tool/module failure, aborting." >&2; exit 1; }
}
run 1000  d1000a
run 1000  d1000b
run 50000 d50000a
run 50000 d50000b

q() { bcftools query -r "chrX:$2-$2" -f '%QUAL\n' "$OUT/repro.$1.vcf.gz" 2>/dev/null | head -1; }

echo "==================================================================="
echo "QUAL at the 7 flip sites — 4 runs"
echo "  same-depth columns EQUAL  => deterministic (real -d effect)"
echo "  same-depth columns DIFFER => nondeterministic (flips are RNG noise)"
echo "==================================================================="
printf "%-10s %-11s %-11s | %-11s %-11s\n" "POS" "d1000a" "d1000b" "d50000a" "d50000b"
for p in $POS; do
  printf "%-10s %-11s %-11s | %-11s %-11s\n" "$p" \
    "$(q d1000a "$p")" "$(q d1000b "$p")" "$(q d50000a "$p")" "$(q d50000b "$p")"
done

echo
echo "Full site record (VCF cols 1-8, INFO included) d=1000a vs d=50000a —"
echo "shows which field carries the QUAL swing at the identical-count sites:"
for p in $POS; do
  echo "--- chrX:$p ---"
  printf "  d1000 : %s\n" "$(bcftools view -H -r "chrX:$p-$p" "$OUT/repro.d1000a.vcf.gz"  2>/dev/null | cut -f1-8 | head -1)"
  printf "  d50000: %s\n" "$(bcftools view -H -r "chrX:$p-$p" "$OUT/repro.d50000a.vcf.gz" 2>/dev/null | cut -f1-8 | head -1)"
done
echo
echo "done."
