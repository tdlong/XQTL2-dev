#!/usr/bin/env bash
# explain_one_snp.sh — for ONE differing SNP, show WHY the tiled caller differs.
#
# Question: at a SNP where RefAlt(validated) != RefAlt(tiled), which SAMPLES
# disagree, and is the disagreement confined to samples whose TRUE (un-capped)
# depth is > 1000 — i.e. exactly the samples the shared "-d 1000" downsamples?
#
#   - disagreement only at true-depth > 1000  -> it's -d 1000 downsampling;
#     quantify the allele-frequency error (expected < 1%). Benign.
#   - disagreement at a sample with true-depth < 1000 -> NOT downsampling.
#     A real bug. Hypothesis falsified; keep digging.
#
# Everything but the TRUE depth is read from the existing RefAlt tables; the
# only compute is one small un-capped mpileup over a 10 kb window at the SNP.
#
# Run from repo root on the cluster:
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/explain_one_snp.sh [chr] [pos]
# (no pos -> auto-picks the first SNP present in BOTH tables with different counts)
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true

REF_DIR=process/AGE_SY_v3
TIL_DIR=process/AGE_SY_v3_tiled
REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
CHR=${1:-chr2R}
POS=${2:-}

a="$REF_DIR/RefAlt.${CHR}.txt"
b="$TIL_DIR/RefAlt.${CHR}.txt"
for f in "$a" "$b" "$REF" "$BAMS"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

# --- pick a SNP present in BOTH tables but with different counts -------------
if [[ -z "$POS" ]]; then
  POS=$(awk 'NR==FNR{ if(FNR>1) A[$2]=$0; next } { if(FNR>1 && ($2 in A) && A[$2]!=$0) {print $2; exit} }' "$a" "$b")
  [[ -z "$POS" ]] && { echo "No shared-but-differing SNP on $CHR (all diffs are presence-only). Try another chr."; exit 0; }
fi
echo "=================================================================="
echo "SNP under the microscope: ${CHR}:${POS}"
echo "=================================================================="

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# sample names (REF_<name> columns in the header)
head -1 "$a" | tr '\t' '\n' | sed -n 's/^REF_//p' | awk '{printf "%03d %s\n", NR, $0}' > "$TMP/names"

# validated caller counts at POS:  idx ref alt   (body: chr pos ref1 alt1 ref2 alt2 ...)
awk -v p="$POS" '$2==p{ for(k=1;k<=(NF-2)/2;k++) printf "%03d %s %s\n", k, $(2*k+1), $(2*k+2); exit }' "$a" > "$TMP/ref"
# tiled caller counts at POS
awk -v p="$POS" '$2==p{ for(k=1;k<=(NF-2)/2;k++) printf "%03d %s %s\n", k, $(2*k+1), $(2*k+2); exit }' "$b" > "$TMP/til"

# TRUE un-capped counts at POS: one small mpileup, -d 100000 (effectively no cap)
P0=$(( POS>5000 ? POS-5000 : 1 )); P1=$(( POS+5000 ))
echo "running un-capped mpileup over ${CHR}:${P0}-${P1} (-d 100000) ..."
bcftools mpileup -I -d 100000 -a FORMAT/AD -f "$REF" -b "$BAMS" -r "${CHR}:${P0}-${P1}" 2>/dev/null \
  | bcftools query -f '%POS[\t%AD]\n' \
  | awk -v p="$POS" '$1==p{ for(k=2;k<=NF;k++){ n=split($k,ad,","); printf "%03d %s %s\n", (k-1), ad[1], (n>=2?ad[2]:0) } exit }' \
  > "$TMP/truth"
[[ -s "$TMP/truth" ]] || { echo "WARN: no mpileup record at POS (SNP not covered in window?)"; }

# --- join everything by sample idx and judge -------------------------------
awk -v NF_names="$TMP/names" -v NF_ref="$TMP/ref" -v NF_til="$TMP/til" -v NF_tru="$TMP/truth" '
  FILENAME==NF_names { name[$1]=$2; next }
  FILENAME==NF_ref   { rr[$1]=$2; ra[$1]=$3; next }
  FILENAME==NF_til   { tr[$1]=$2; ta[$1]=$3; next }
  FILENAME==NF_tru   { ur[$1]=$2; ua[$1]=$3; next }
  END{
    printf "%-16s %10s %10s %10s   %6s %6s   %s\n","sample","validated","tiled","TRUEdepth","vFreq","tFreq","note"
    ndis=0; maxd=0; badlow=0
    ns=0; for(k in name) ns++
    for(i=1;i<=ns;i++){
      k = sprintf("%03d", i)
      vr=rr[k]+0; va=ra[k]+0; xr=tr[k]+0; xa=ta[k]+0; dr=ur[k]+0; da=ua[k]+0; dp=dr+da
      disagree = (vr!=xr || va!=xa)
      vf = (vr+va>0)? vr/(vr+va):0
      xf = (xr+xa>0)? xr/(xr+xa):0
      note=""
      if(disagree){ ndis++; d=vf-xf; if(d<0)d=-d; if(d>maxd)maxd=d
                    note=sprintf("DISAGREE dFreq=%+.4f", vf-xf)
                    if(dp>0 && dp<=1000){ note=note "  <<< TRUE DEPTH <=1000 : NOT downsampling!"; badlow++ }
                    else note=note "  (true depth "dp" > 1000: consistent w/ -d cap)" }
      printf "%-16s  %4d/%-4d  %4d/%-4d  %9d   %.3f  %.3f   %s\n", name[k], vr,va, xr,xa, dp, vf,xf, note
    }
    print  "------------------------------------------------------------------"
    printf "samples disagreeing: %d   max |dFreq|: %.4f\n", ndis, maxd
    if(badlow>0) printf "VERDICT: %d disagreeing sample(s) had TRUE depth <=1000 -> NOT downsampling. Real bug.\n", badlow
    else         printf "VERDICT: every disagreement is at a sample with TRUE depth >1000 -> consistent with -d 1000 downsampling (freq error <= %.2f%%).\n", maxd*100
  }
' "$TMP/names" "$TMP/ref" "$TMP/til" "$TMP/truth"
