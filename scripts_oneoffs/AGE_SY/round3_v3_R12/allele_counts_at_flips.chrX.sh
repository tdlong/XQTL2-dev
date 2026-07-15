#!/usr/bin/env bash
# allele_counts_at_flips.chrX.sh — dump the ENTIRE per-sample table for the 7
# flip sites. No summaries. For each SNP, every one of the 72 samples, its depth
# and its allele read counts, at d=1000 and d=50000, side by side.
#
# Reading a row: DP = reads at that site in that sample; AD = comma-separated
# reads per allele in the order shown in the SNP header (ref first). If a
# sample's DP is the SAME at d=1000 and d=50000 it was never capped; if d=1000
# DP is 1000 and d=50000 is larger, that sample was subsampled.
#
# Run from repo root on the cluster (a few minutes):
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/allele_counts_at_flips.chrX.sh
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true

REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
POSITIONS="16025382 16030427 16030958 16031180 16031851 16045058 16363526"
OUT=process/flip_detail; mkdir -p "$OUT"
REG="$OUT/regions.txt"; : > "$REG"
for p in $POSITIONS; do printf 'chrX\t%s\n' "$p" >> "$REG"; done

# sample names, in the order the columns appear
bcftools mpileup -I -d 50 -r "chrX:${POSITIONS%% *}-${POSITIONS%% *}" -f "$REF" -b "$BAMS" 2>/dev/null \
  | bcftools query -l > "$OUT/samples.txt"

# raw capture: POS REF ALT QUAL then, per sample, "DP:AD"
capture() {  # $1=depth  $2=outfile
  echo "  mpileup -d $1 at the 7 sites ..." >&2
  bcftools mpileup -I -d "$1" -R "$REG" -a FORMAT/AD,FORMAT/DP -f "$REF" -b "$BAMS" 2>/dev/null \
    | bcftools call -mv 2>/dev/null \
    | bcftools query -f '%POS\t%REF\t%ALT\t%QUAL[\t%DP:%AD]\n' > "$2"
}
capture 1000  "$OUT/d1000.tsv"
capture 50000 "$OUT/d50000.tsv"
echo "raw saved (no re-run needed): $OUT/d1000.tsv  $OUT/d50000.tsv  $OUT/samples.txt"
echo

awk -v POS="$POSITIONS" -v SN="$OUT/samples.txt" -v C="$OUT/d1000.tsv" -v U="$OUT/d50000.tsv" '
  function verdict(ref,alt,q,   nf,aa,snp,j){ if(alt=="")return "no-call"
    nf=split(alt,aa,","); snp=(length(ref)==1); for(j=1;j<=nf;j++) if(length(aa[j])!=1)snp=0
    return (nf==1 && snp && q+0>59)?"kept":"THROWN-OUT" }
  BEGIN{
    while((getline l < SN)>0) s[++ns]=l
    while((getline l < C)>0){ n=split(l,f,"\t"); p=f[1]; rc[p]=f[2]; ac[p]=f[3]; qc[p]=f[4]; inC[p]=1
                              for(i=5;i<=n;i++) cc[p,i-4]=f[i] }
    while((getline l < U)>0){ n=split(l,f,"\t"); p=f[1]; ru[p]=f[2]; au[p]=f[3]; qu[p]=f[4]; inU[p]=1
                              for(i=5;i<=n;i++) cu[p,i-4]=f[i] }
    np=split(POS,plist," ")
    for(x=1;x<=np;x++){ p=plist[x]
      print  "############################################################################"
      printf "chrX:%s   REF=%s\n", p, (inC[p]?rc[p]:(inU[p]?ru[p]:"?"))
      printf "  d=1000 : ALT=%-8s QUAL=%-8s -> %s\n", (inC[p]?ac[p]:"-"), (inC[p]?qc[p]:"-"), (inC[p]?verdict(rc[p],ac[p],qc[p]):"no-call")
      printf "  d=50000: ALT=%-8s QUAL=%-8s -> %s\n", (inU[p]?au[p]:"-"), (inU[p]?qu[p]:"-"), (inU[p]?verdict(ru[p],au[p],qu[p]):"no-call")
      printf "  %-20s | %-18s | %-18s\n", "sample", "d1000 DP:AD", "d50000 DP:AD"
      for(k=1;k<=ns;k++)
        printf "  %-20s | %-18s | %-18s\n", s[k], (inC[p]?cc[p,k]:"-"), (inU[p]?cu[p,k]:"-")
    }
  }'
