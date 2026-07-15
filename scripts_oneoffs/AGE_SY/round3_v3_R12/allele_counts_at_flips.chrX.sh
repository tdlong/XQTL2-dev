#!/usr/bin/env bash
# allele_counts_at_flips.chrX.sh — the ACTUAL numbers behind the 7 flipped sites.
#
# Fetches ONLY the 7 positions (-r, index jump) and calls them at d=1000 and
# d=50000. For each, prints: REF>ALT, QUAL, kept/thrown, and total reads
# supporting ref / alt1 / alt2 ... summed over all 72 samples. Fast (minutes).
#
# These verdicts should match the earlier 7-flip table; if any differ, the call
# is window-sensitive and that itself is worth knowing.
#
# Run from repo root on the cluster (a few minutes, no sbatch needed):
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/allele_counts_at_flips.chrX.sh
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true

REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
POSITIONS="16025382 16030427 16030958 16031180 16031851 16045058 16363526"

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
for p in $POSITIONS; do printf 'chrX\t%s\n' "$p"; done > "$TMP/regions"

# per position -> "REF>ALT  QUAL=..  kept|THROWN-OUT  ref=N alt1=N alt2=N" (AD summed over samples)
dump() {  # $1 = depth   $2 = outfile
  echo "  calling the 7 sites at -d $1 ..." >&2
  bcftools mpileup -I -d "$1" -R "$TMP/regions" -a FORMAT/AD -f "$REF" -b "$BAMS" 2>/dev/null \
    | bcftools call -mv 2>/dev/null \
    | bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t[%AD ]\n' \
    | awk '{
        delete tot; nall=0
        for(i=5;i<=NF;i++){ n=split($i,a,","); if(n>nall)nall=n
                            for(j=1;j<=n;j++) if(a[j]!=".") tot[j]+=a[j] }
        cs="ref="tot[1]; for(j=2;j<=nall;j++) cs=cs"  alt"(j-1)"="tot[j]
        nalt=split($3,aa,","); issnp=(length($2)==1)
        for(j=1;j<=nalt;j++) if(length(aa[j])!=1) issnp=0
        verdict=(nalt==1 && issnp && $4>59)?"kept":"THROWN-OUT"
        printf "%s\t%s>%s\tQUAL=%s\t%-10s\t%s\n", $1,$2,$3,$4,verdict,cs
      }' > "$2"
}

dump 1000  "$TMP/d1000"
dump 50000 "$TMP/d50000"

echo
echo "==================================================================="
echo "The 7 flipped sites — REF>ALT, QUAL, verdict, total reads per allele"
echo "(read counts summed over all 72 samples; 'kept' = passes biallelic + QUAL>59)"
echo "==================================================================="
for p in $POSITIONS; do
  a=$(awk -v p="$p" '$1==p{sub(/^[0-9]+\t/,""); print}' "$TMP/d1000")
  b=$(awk -v p="$p" '$1==p{sub(/^[0-9]+\t/,""); print}' "$TMP/d50000")
  echo "chrX:$p"
  printf "   d=1000    %s\n" "${a:-not called here (no variant)}"
  printf "   d=50000   %s\n" "${b:-not called here (no variant)}"
  echo
done
