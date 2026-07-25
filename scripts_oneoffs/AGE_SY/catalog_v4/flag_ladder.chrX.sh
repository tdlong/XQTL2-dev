#!/bin/bash
#SBATCH --job-name=flag_ladder
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=03:00:00
#SBATCH -o logs/AGE_SY/flag_ladder_%j.out
###############################################################################
# flag_ladder.chrX.sh — PER-SITE re-count of one sample over the chrX catalog
# under the full q/Q/BAQ flag matrix, SAME caller (call -m -C alleles), so only
# the mpileup flag varies. Writes one wide table (ref/alt per setting per site)
# so we can compare CALLS site-by-site, not just summed totals.
#
# 2x2 of the two mpileup filters, all BAQ-off (that is the catalog method — BAQ is
# NOT toggled here; the BAQ-on reference is the real v3 table, compared separately):
#   v4    : -B -q20 -Q20     -q off : -B -q0 -Q20     -Q off : -B -q20 -Q13     pivot : -B -q0 -Q13
#
# OUTPUT: process/AGE_SY_v4/flag_ladder.<sample>.chrX.tsv.gz
#   POS  ref.v4 alt.v4  ref.qoff alt.qoff  ref.Qoff alt.Qoff  ref.pivot alt.pivot
#
#   sbatch scripts_oneoffs/AGE_SY/catalog_v4/flag_ladder.chrX.sh [sample_bam]
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true
REF=pipeline/ref/dm6.fa
CAT=process/AGE_SY_v4/Catalog/catalog.tsv.gz
BAM=${1:-data/bam/AGE_SY/Con_R5_F.bam}
SAMP=$(basename "$BAM" .bam)
OUT=process/AGE_SY_v4/flag_ladder.${SAMP}.chrX.tsv.gz
for f in "$REF" "$CAT" "$BAM"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
echo "sample $SAMP  chrX  -> $OUT"

# per-setting: emit "ref<TAB>alt" per catalog chrX site (sites in catalog order, aligned across runs)
call() {  # $1=outfile  $2..=mpileup flags
  local out="$1"; shift
  bcftools mpileup "$@" -r chrX -T "$CAT" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
    | bcftools call -m -C alleles -T "$CAT" 2>/dev/null \
    | bcftools query -f '%POS\t[%AD]\n' \
    | awk -F'\t' '{split($2,a,","); print (a[1]==""||a[1]=="."?0:a[1])"\t"(a[2]==""||a[2]=="."?0:a[2])}' > "$out"
}
bcftools mpileup -B -q0 -Q13 -r chrX -T "$CAT" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
  | bcftools call -m -C alleles -T "$CAT" 2>/dev/null | bcftools query -f '%POS\n' > "$TMP/pos"
call "$TMP/v4"     -B -q 20 -Q 20
call "$TMP/qoff"   -B -q 0  -Q 20
call "$TMP/Qoff"   -B -q 20 -Q 13
call "$TMP/pivot"  -B -q 0  -Q 13

{ echo -e "POS\tref.v4\talt.v4\tref.qoff\talt.qoff\tref.Qoff\talt.Qoff\tref.pivot\talt.pivot"
  paste "$TMP/pos" "$TMP/v4" "$TMP/qoff" "$TMP/Qoff" "$TMP/pivot"
} | bgzip > "$OUT"
echo "rows: $(zcat "$OUT" | wc -l)   wrote $OUT"
echo "summed ALT% per setting (quick check):"
zcat "$OUT" | awk -F'\t' 'NR>1{v4r+=$2;v4a+=$3;qr+=$4;qa+=$5;Qr+=$6;Qa+=$7;pr+=$8;pa+=$9}
  END{printf "  v4(-q20-Q20)=%.3f  -qoff(-q0-Q20)=%.3f  -Qoff(-q20-Q13)=%.3f  pivot(-q0-Q13)=%.3f  (ALT%%)\n",
      100*v4a/(v4r+v4a),100*qa/(qr+qa),100*Qa/(Qr+Qa),100*pa/(pr+pa)}'
