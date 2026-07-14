#!/usr/bin/env bash
# profile_diff_snps.sh — true (un-capped) coverage at the SNPs where the two
# callers disagree, vs a control set where they agree.
#
# Hypothesis: the disagreements are the -d 1000 downsampling, so they should sit
# at high-coverage sites (some sample's true depth > 1000). Agreeing SNPs should
# not. One mpileup over each position set (index-jump via -R), -d 100000 so
# nothing is capped.
#
#   differing SNPs mostly high-cov + control mostly low-cov -> it's the cap.
#   differing SNPs NOT enriched for high-cov               -> a real bug.
#
# Run from repo root on the cluster:
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/profile_diff_snps.sh [ "chrX chr2L ..." ]
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true

REF_DIR=process/AGE_SY_v3
TIL_DIR=process/AGE_SY_v3_tiled
REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
CHRS=${1:-"chrX chr2L chr2R chr3L chr3R"}
CAP=1000

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
: > "$TMP/diff.rgn"; : > "$TMP/ctrl.rgn"

for chr in $CHRS; do
  a="$REF_DIR/RefAlt.$chr.txt"; b="$TIL_DIR/RefAlt.$chr.txt"
  [[ -f "$a" && -f "$b" ]] || { echo "skip $chr (missing table)"; continue; }
  # positions where the callers' lines differ (present-in-one counts as differ)
  awk -v chr="$chr" 'NR==FNR{ if(FNR>1)A[$2]=$0; next }
                     FNR>1{ B[$2]=1; if(!($2 in A) || A[$2]!=$0) print chr,$2 }
                     END{ for(p in A) if(!(p in B)) print chr,p }' OFS='\t' "$a" "$b" >> "$TMP/diff.rgn"
  # control: positions where the callers' lines are identical (random 200)
  awk -v chr="$chr" 'NR==FNR{ if(FNR>1)A[$2]=$0; next }
                     FNR>1 && ($2 in A) && A[$2]==$0 { print chr,$2 }' OFS='\t' "$a" "$b" \
    | shuf | head -200 >> "$TMP/ctrl.rgn"
done
sort -k1,1 -k2,2n -o "$TMP/diff.rgn" "$TMP/diff.rgn"
sort -k1,1 -k2,2n -o "$TMP/ctrl.rgn" "$TMP/ctrl.rgn"

# per SNP: max per-sample true depth, and how many samples exceed the cap
depth_profile () {
  bcftools mpileup -I -d 100000 -a FORMAT/DP -f "$REF" -b "$BAMS" -R "$1" 2>/dev/null \
   | bcftools query -f '%CHROM\t%POS[\t%DP]\n' \
   | awk -v cap="$CAP" '{ mx=0; nov=0; for(i=3;i<=NF;i++){ d=$i+0; if(d>mx)mx=d; if(d>cap)nov++ }
                          n++; if(mx>cap)hi++; sum+=mx; if(mx>mmax)mmax=mx; sumov+=nov }
       END{ if(n==0){print "  (no SNPs)"; next}
            printf "  %d SNPs | %d (%.0f%%) have a sample over %d cov | meanMaxDepth=%.0f | maxMaxDepth=%d | avg #capped-samples/SNP=%.1f\n",
                   n, hi, 100*hi/n, cap, sum/n, mmax, sumov/n }'
}

echo "=================================================================="
echo "differing SNPs : $(wc -l < "$TMP/diff.rgn")     control (agreeing) SNPs: $(wc -l < "$TMP/ctrl.rgn")"
echo "=================================================================="
echo "TRUE (un-capped) coverage at DIFFERING SNPs:"
depth_profile "$TMP/diff.rgn"
echo "TRUE (un-capped) coverage at CONTROL (agreeing) SNPs:"
depth_profile "$TMP/ctrl.rgn"
echo "------------------------------------------------------------------"
echo "If differing SNPs are ~all high-cov and controls are not -> -d 1000 cap."
echo "If differing SNPs are NOT enriched for high-cov          -> a real bug."
