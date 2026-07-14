#!/usr/bin/env bash
# characterize_diff_snps.sh — plain description of how the two callers' SNP
# tables differ. No mpileup, just arithmetic on the existing RefAlt tables.
#
# OLD = validated caller (process/AGE_SY_v3)
# NEW = tiled caller     (process/AGE_SY_v3_tiled)
#
# Answers:
#   1. When BOTH callers report a position, are the 72-sample counts identical?
#   2. How many positions are reported only by OLD? only by NEW?
#   3. For each such position: total reads over all samples, and the
#      alternate-allele frequency over all samples (sum ALT / sum REF+ALT).
#
# Run from repo root on the cluster:
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/characterize_diff_snps.sh [chr ...]
set -uo pipefail
REF_DIR=process/AGE_SY_v3
TIL_DIR=process/AGE_SY_v3_tiled
CHRS=${1:-"chrX chr2L chr2R chr3L chr3R"}

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
: > "$TMP/difflist"
ts=0; tc=0; tr=0; tt=0

echo "How the two tables compare (one row per SNP position):"
printf "  %-7s %12s %12s %10s %10s\n" "chr" "both_same" "both_DIFFER" "OLD_only" "NEW_only"
for chr in $CHRS; do
  a="$REF_DIR/RefAlt.$chr.txt"; b="$TIL_DIR/RefAlt.$chr.txt"
  [[ -f "$a" && -f "$b" ]] || { echo "  skip $chr (missing table)"; continue; }
  read s c r t < <(awk -v chr="$chr" -v out="$TMP/difflist" '
    function tf(line,   f,nf,k,sr,sa){ nf=split(line,f," "); sr=0; sa=0
      for(k=1;k<=(nf-2)/2;k++){ sr+=f[2*k+1]; sa+=f[2*k+2] }
      return (sr+sa)"\t"(sr+sa>0 ? sa/(sr+sa) : 0) }
    NR==FNR{ if(FNR>1) A[$2]=$0; next }
    FNR>1 { B[$2]=$0 }
    END{ same=0;chg=0;ro=0;to=0
      for(p in A){ if(p in B){ if(A[p]==B[p]) same++
                               else { chg++; printf "%s\t%s\t%s\t%s\n",chr,p,"BOTH-DIFFER",tf(A[p]) >> out } }
                   else { ro++; printf "%s\t%s\t%s\t%s\n",chr,p,"OLD-only",tf(A[p]) >> out } }
      for(p in B){ if(!(p in A)){ to++; printf "%s\t%s\t%s\t%s\n",chr,p,"NEW-only",tf(B[p]) >> out } }
      print same, chg, ro, to }' "$a" "$b")
  printf "  %-7s %12d %12d %10d %10d\n" "$chr" "$s" "$c" "$r" "$t"
  ts=$((ts+s)); tc=$((tc+c)); tr=$((tr+r)); tt=$((tt+t))
done
printf "  %-7s %12d %12d %10d %10d\n" "TOTAL" "$ts" "$tc" "$tr" "$tt"

echo
echo "FACT: positions reported by BOTH callers but with DIFFERENT counts = $tc"
[ "$tc" -eq 0 ] && echo "      -> zero. Whenever both callers report a position, all 72-sample counts are IDENTICAL."
echo
echo "Every differing SNP (present in only one caller), sorted by alt-allele frequency:"
echo "  chr  pos  caller  total_reads(all samples)  alt_freq(all samples)"
sort -t$'\t' -k5,5g "$TMP/difflist" \
  | awk -F'\t' '{printf "  %-6s %-10s %-9s reads=%-9s altfreq=%.4f\n",$1,$2,$3,$4,$5}'
