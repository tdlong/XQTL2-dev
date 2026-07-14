#!/usr/bin/env bash
# counts_at_diff_snps.sh — read the ACTUAL per-sample counts at the differing
# SNPs straight from the RefAlt tables (no mpileup) and test the -d 1000 story.
#
# Every difference is present-in-one-caller (a SNP one caller emits, the other
# doesn't). The emitting caller stored its per-sample ref/alt counts. If the
# cause is the -d 1000 cap, some sample at these SNPs must be pegged: ref+alt
# ~ 1000. Agreeing SNPs (control) should not be pegged.
#
#   PREDICTION (cap):  differing SNPs -> a sample with ref+alt ~1000; controls not.
#   FALSIFIED (bug):   differing SNPs have all samples well under 1000.
#
# Run from repo root on the cluster:
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/counts_at_diff_snps.sh [chr]
set -uo pipefail

REF_DIR=process/AGE_SY_v3
TIL_DIR=process/AGE_SY_v3_tiled
CHRS=${1:-"chrX chr2L chr2R chr3L chr3R"}
PEG=990          # ref+alt at/above this = effectively capped at -d 1000

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
: > "$TMP/diff.tsv"; : > "$TMP/ctrl.tsv"     # rows: chr pos src maxTotal nPegged

for chr in $CHRS; do
  a="$REF_DIR/RefAlt.$chr.txt"; b="$TIL_DIR/RefAlt.$chr.txt"
  [[ -f "$a" && -f "$b" ]] || { echo "skip $chr (missing table)"; continue; }
  awk -v chr="$chr" -v peg="$PEG" -v d="$TMP/diff.tsv" -v c="$TMP/ctrl.full" '
    function stats(line, src, out,   f,nf,k,t,mx,np){
      nf=split(line,f," "); mx=0; np=0
      for(k=1;k<=(nf-2)/2;k++){ t=f[2*k+1]+f[2*k+2]; if(t>mx)mx=t; if(t>=peg)np++ }
      printf "%s\t%s\t%s\t%d\t%d\n", chr, f[2], src, mx, np >> out
    }
    NR==FNR{ if(FNR>1) A[$2]=$0; next }
    FNR>1  { B[$2]=$0 }
    END{
      for(p in A){ if(!(p in B)) stats(A[p],"REFonly",d)
                   else if(A[p]!=B[p]) stats(A[p],"CHANGED",d) }
      for(p in B){ if(!(p in A)) stats(B[p],"TILEDonly",d) }
      for(p in A){ if((p in B) && A[p]==B[p]) stats(A[p],"agree",c) }
    }' "$a" "$b"
done
# control: random 300 agreeing SNPs
shuf "$TMP/ctrl.full" 2>/dev/null | head -300 > "$TMP/ctrl.tsv" || true

summ () {   # $1 = tsv  $2 = label
  awk -v lab="$2" -v peg="$PEG" '
    { n++; if($5>0)withpeg++; sum+=$4; if($4>mx)mx=$4 }
    END{ if(n==0){ printf "  %-24s (none)\n", lab; }
         else printf "  %-24s %5d SNPs | %5.1f%% have a sample with ref+alt>=%d | meanMaxTotal=%.0f | maxMaxTotal=%d\n",
                     lab, n, 100*withpeg/n, peg, sum/n, mx }' "$1"
}

echo "=================================================================="
echo "Example differing SNPs (chr pos source maxTotal nSamplesPegged):"
sort -k4,4nr "$TMP/diff.tsv" | head -8 | awk '{printf "  %-6s %-10s %-10s maxTotal=%-5s peggedSamples=%s\n",$1,$2,$3,$4,$5}'
echo "------------------------------------------------------------------"
echo "Per-sample count peg (ref+alt), from the existing tables:"
summ "$TMP/diff.tsv" "DIFFERING SNPs:"
summ "$TMP/ctrl.tsv" "CONTROL (agreeing) SNPs:"
echo "------------------------------------------------------------------"
echo "Differing pegged ~1000 AND controls not -> the -d 1000 cap."
echo "Differing NOT pegged                    -> cap is not the cause; real bug."
