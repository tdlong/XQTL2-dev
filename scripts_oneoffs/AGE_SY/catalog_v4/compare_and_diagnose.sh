#!/bin/bash
#SBATCH --job-name=cat_v4_compare
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:40:00
#SBATCH -o logs/AGE_SY/compare_v3_v4_%j.out
###############################################################################
# compare_and_diagnose.sh — verify the founder-catalog v4 callset vs v3.
#
# New XQTL2 layout (d40e279): catalog in process/AGE_SY_v4/Catalog, calls in
# process/AGE_SY_v4/Calls. Answers the three things XQTL2-dev #1 asks for:
#   0. positions unique (the #14 symptom is gone)  -> expect 0
#   1. per-chr SNP counts, v3 vs v4
#   2. catalog.stats.txt: per-rule SNP tally + kept size vs v3
#   3. compare_refalt_calls.R: counts agree at shared SNPs? (clean, no dcast warnings)
#
# Submit from repo root (usually auto-chained after call_samples by run_catalog_v4.sh):
#   sbatch scripts_oneoffs/AGE_SY/catalog_v4/compare_and_diagnose.sh
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true

V3=process/AGE_SY_v3                 # validated QUAL caller
V4=process/AGE_SY_v4/Calls           # candidate founder-catalog calls
CAT=process/AGE_SY_v4/Catalog        # the catalog
CHRS="chrX chr2L chr2R chr3L chr3R"

echo "=================================================================="
echo "0) DUPLICATION CHECK (XQTL2 #14): duplicate (CHROM,POS) per chr — expect 0"
echo "=================================================================="
for c in $CHRS; do
  f="$V4/RefAlt.$c.txt"
  [[ -e "$f" ]] || { printf "  %-7s (missing %s)\n" "$c" "$f"; continue; }
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

echo
echo "=================================================================="
echo "2) catalog.stats.txt — per-rule SNP tally (candidates, dropped, kept)"
echo "=================================================================="
if [[ -e "$CAT/catalog.stats.txt" ]]; then cat "$CAT/catalog.stats.txt"; else echo "  (missing $CAT/catalog.stats.txt)"; fi

echo
echo "=================================================================="
echo "3) compare_refalt_calls.R: v3 vs v4 — counts agree at shared SNPs?"
echo "   (should run clean now — no dcast duplicate-key warnings)"
echo "   NOTE: if v4 has fewer samples than v3 (phase 1), read SITE overlap,"
echo "   not per-sample count agreement (that needs the full 72 in v4)."
echo "=================================================================="
module load R/4.2.2 2>/dev/null || true
Rscript pipeline/scripts/compare_refalt_calls.R "$V3" "$V4" || \
  echo "  (compare_refalt_calls.R exited nonzero)"
echo
echo "done."
