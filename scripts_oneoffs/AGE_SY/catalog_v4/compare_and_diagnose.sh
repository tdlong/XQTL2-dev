#!/bin/bash
#SBATCH --job-name=cat_v4_compare
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:40:00
#SBATCH -o logs/AGE_SY/compare_v3_v4_%j.out
###############################################################################
# compare_and_diagnose.sh — one job: evaluate the v4 founder-catalog caller
# against v3 (validated QUAL caller), AND diagnose the chr2L shortfall.
#
# Phase-1 v4 has 36 samples (R1-R6); v3 has 72. So the per-SAMPLE count
# comparison is not apples-to-apples yet (that needs phase 2). Here we read the
# parts that ARE meaningful now: per-chr SNP counts, SNP-site overlap, and the
# chr2L cause (catalog site count is set by the founders only, sample-independent).
#
# Submit from repo root:
#   sbatch scripts_oneoffs/AGE_SY/catalog_v4/compare_and_diagnose.sh
#   -> process/AGE_SY_v4/compare_v3_v4_<jobid>.out
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true   # NO samtools (issue #13 clash)

V3=process/AGE_SY_v3
V4=process/AGE_SY_v4
CHRS="chrX chr2L chr2R chr3L chr3R"

echo "=================================================================="
echo "0) DUPLICATION CHECK (confirms XQTL2 #14 fix): duplicate (CHROM,POS) per chr"
echo "   expect 0 everywhere; pre-fix chrX had 66,793."
echo "=================================================================="
for c in $CHRS; do
  f="$V4/RefAlt.$c.txt"
  [[ -e "$f" ]] || { printf "  %-7s (missing)\n" "$c"; continue; }
  d=$(tail -n +2 "$f" | cut -f1,2 | sort | uniq -d | wc -l | tr -d ' ')
  printf "  %-7s duplicate positions: %s\n" "$c" "$d"
done

echo
echo "=================================================================="
echo "1) per-chr SNP counts: v3 (QUAL) vs v4 (founder-catalog)"
echo "=================================================================="
printf "  %-7s %12s %12s %8s\n" chr v3 v4 "v4/v3"
for c in $CHRS; do
  a=$(wc -l < "$V3/RefAlt.$c.txt" 2>/dev/null || echo 0)
  b=$(wc -l < "$V4/RefAlt.$c.txt" 2>/dev/null || echo 0)
  r=$(awk -v a="$a" -v b="$b" 'BEGIN{printf (a>0)?"%.2f":"NA", b/a}')
  printf "  %-7s %12s %12s %8s\n" "$c" "$a" "$b" "$r"
done
echo "  (chr2L is the arm under scrutiny — expect it to look ~8x low in v4)"

echo
echo "=================================================================="
echo "2) chr2L cause: per-founder mean depth, chr2L vs chr2R"
echo "   (catalog needs EVERY founder covered; a founder that drops out on"
echo "    chr2L would gate out the whole arm. Reads the founder BCFs on disk.)"
echo "=================================================================="
for c in chr2L chr2R; do
  bcf="$V4/founders.calls.$c.bcf"
  [[ -e "$bcf" ]] || { echo "  missing $bcf"; continue; }
  echo "  --- $c (mean DP per founder over first 50k sites) ---"
  # header sample order
  bcftools query -l "$bcf" | tr '\n' '\t' | sed 's/^/    samples: /'; echo
  bcftools query -f '[%DP\t]\n' "$bcf" 2>/dev/null | head -50000 \
    | awk -F'\t' '{ for(i=1;i<=NF;i++) if($i!=""){ s[i]+=($i=="."?0:$i); n[i]++ } }
                  END{ printf "    meanDP :"; for(i=1;i<=length(s);i++) printf " %7.1f", s[i]/n[i]; printf "\n" }'
done
echo "  -> a founder whose meanDP collapses on chr2L (vs chr2R) is the culprit;"
echo "     if ALL founders drop, it is a chr2L reference/mapping issue, not one BAM."

echo
echo "=================================================================="
echo "3) compare_refalt_calls.R: v3 vs v4 (SNP-site overlap)"
echo "   NOTE phase-1 sample mismatch (36 vs 72) — read the SITE overlap,"
echo "   not the per-sample count agreement (that needs phase 2)."
echo "=================================================================="
module load R/4.2.2 2>/dev/null || true
Rscript pipeline/scripts/compare_refalt_calls.R "$V3" "$V4" || \
  echo "  (compare_refalt_calls.R exited nonzero — likely the 36-vs-72 sample mismatch; phase 2 fixes this)"
echo
echo "done."
