#!/bin/bash
#SBATCH --job-name=baq_caller
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=03:00:00
#SBATCH -o logs/AGE_SY/baq_caller_%j.out
###############################################################################
# baq_caller.chrX.sh — the flag ladder showed -q/-Q do NOT change ALT. So the
# v3-vs-v4 ALT gap is BAQ or the caller constraint. Test exactly those, 2x2, with
# the CORRECT --max-depth 2000 (my earlier ladder defaulted to 250 and undercounted).
# All: -q0 -Q13 --max-depth 2000 -T catalog, one sample, chrX, per-site table.
#
#   catalog : -B      call -m -C alleles   (BAQ off, constrained = real v4 recipe)
#   +BAQ    :         call -m -C alleles   (BAQ on,  constrained)   -> BAQ effect
#   +free   : -B      call -m               (BAQ off, unconstrained)-> constraint effect
#   v3ish   :         call -m               (BAQ on,  unconstrained ~ v3)
#
#   sbatch scripts_oneoffs/AGE_SY/catalog_v4/baq_caller.chrX.sh [sample_bam]
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true
REF=pipeline/ref/dm6.fa
CAT=process/AGE_SY_v4/Catalog/catalog.tsv.gz
BAM=${1:-data/bam/AGE_SY/Con_R5_F.bam}
SAMP=$(basename "$BAM" .bam)
OUT=process/AGE_SY_v4/baq_caller.${SAMP}.chrX.tsv.gz
for f in "$REF" "$CAT" "$BAM"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
echo "sample $SAMP  chrX  --max-depth 2000  -q0 -Q13 -> $OUT"

# ref = AD[1]; alt = (sum of AD) - ref  (robust to multiallelic in the unconstrained calls)
run() {  # $1=out  $2=BAQ("-B" or "")  $3=constrain("-C alleles" or "")
  bcftools mpileup $2 -q 0 -Q 13 --max-depth 2000 -r chrX -T "$CAT" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
    | bcftools call -m $3 -T "$CAT" 2>/dev/null \
    | bcftools query -f '%POS\t[%AD]\n' \
    | awk -F'\t' '{n=split($2,a,","); r=(a[1]=="."?0:a[1]); t=0; for(i=1;i<=n;i++) t+=(a[i]=="."?0:a[i]); print r"\t"(t-r)}' > "$1"
}
bcftools mpileup -B -q0 -Q13 --max-depth 2000 -r chrX -T "$CAT" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
  | bcftools call -m -C alleles -T "$CAT" 2>/dev/null | bcftools query -f '%POS\n' > "$TMP/pos"
run "$TMP/catalog" "-B" "-C alleles"
run "$TMP/baqon"   ""   "-C alleles"
run "$TMP/free"    "-B" ""
run "$TMP/v3ish"   ""   ""

{ echo -e "POS\tref.catalog\talt.catalog\tref.baqon\talt.baqon\tref.free\talt.free\tref.v3ish\talt.v3ish"
  paste "$TMP/pos" "$TMP/catalog" "$TMP/baqon" "$TMP/free" "$TMP/v3ish"
} | bgzip > "$OUT"
echo "rows: $(zcat "$OUT" | wc -l)   wrote $OUT"
zcat "$OUT" | awk -F'\t' 'NR>1{cr+=$2;ca+=$3;br+=$4;ba+=$5;fr+=$6;fa+=$7;vr+=$8;va+=$9}
  END{printf "ALT%%  catalog(-B,-C)=%.3f  +BAQ(-C)=%.3f  +free(-B,-m)=%.3f  v3ish(-m)=%.3f\n",
      100*ca/(cr+ca),100*ba/(br+ba),100*fa/(fr+fa),100*va/(vr+va)}'
echo "read: catalog->+BAQ = BAQ effect; catalog->+free = the -C alleles constraint; v3ish should ~ real v3"
