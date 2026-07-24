#!/bin/bash
#SBATCH --job-name=cat_v4_catcheck
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:15:00
#SBATCH -o logs/AGE_SY/catalog_check_%j.out
###############################################################################
# catalog_check.sh — sanity-check the built catalog BEFORE calling samples.
# Reports the per-rule tally and per-chromosome SNP counts (is chr2L recovered?),
# so we don't call 44 samples against a bad catalog. Read via cluster_sync.
#
#   sbatch scripts_oneoffs/AGE_SY/catalog_v4/catalog_check.sh
###############################################################################
set -uo pipefail
CAT=process/AGE_SY_v4/Catalog
V3=process/AGE_SY_v3

echo "=== catalog.stats.txt (per-rule SNP tally) ==="
[[ -e "$CAT/catalog.stats.txt" ]] && cat "$CAT/catalog.stats.txt" || echo "  (missing $CAT/catalog.stats.txt)"

echo
echo "=== format check: first 3 lines of catalog.tsv.gz (want CHROM POS REF,ALT) ==="
zcat "$CAT/catalog.tsv.gz" 2>/dev/null | head -3

echo
echo "=== catalog SNPs per chromosome vs v3 RefAlt (was chr2L: 15,208 broken) ==="
TMP=$(mktemp)
zcat "$CAT/catalog.tsv.gz" 2>/dev/null | cut -f1 | sort | uniq -c | awk '{print $2"\t"$1}' > "$TMP"
printf "  %-7s %12s %12s %8s\n" chr catalog v3 "cat/v3"
for c in chrX chr2L chr2R chr3L chr3R; do
  n=$(awk -v c="$c" '$1==c{print $2}' "$TMP"); n=${n:-0}
  v=$(wc -l < "$V3/RefAlt.$c.txt" 2>/dev/null || echo 0)
  r=$(awk -v n="$n" -v v="$v" 'BEGIN{printf (v>0)?"%.2f":"NA", n/v}')
  printf "  %-7s %12s %12s %8s\n" "$c" "$n" "$v" "$r"
done
rm -f "$TMP"
echo
echo "total catalog SNPs: $(zcat "$CAT/catalog.tsv.gz" 2>/dev/null | wc -l)"
echo "done."
