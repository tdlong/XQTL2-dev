#!/bin/bash
#SBATCH --job-name=allele_counts
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=6:00:00
#SBATCH -o allele_counts_at_flips.chrX.%j.out
###############################################################################
# allele_counts_at_flips.chrX.sh — the ACTUAL numbers behind the 7 flipped sites.
#
# Re-runs the caller over the same window at d=1000 and d=50000 and, for each of
# the 7 positions where the call changed, prints:
#   REF base, ALT allele(s), QUAL, kept/thrown, and the total reads supporting
#   ref / alt1 / alt2 ... summed over all 72 samples.
# This is what explains WHY each site flips. No interpretation, just the numbers.
#
# Submit from repo root on the cluster (~2h; two full-window mpileups):
#   sbatch scripts_oneoffs/AGE_SY/round3_v3_R12/allele_counts_at_flips.chrX.sh
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true

REGION=chrX:15900000-16400000
REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
POSITIONS="16025382 16030427 16030958 16031180 16031851 16045058 16363526"

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# per position -> "REF  ALT  QUAL  kept|THROWN  ref=N alt1=N alt2=N ..." (AD summed over samples)
dump() {  # $1 = depth   $2 = outfile
  echo "  running mpileup -d $1 over $REGION ..." >&2
  bcftools mpileup -I -d "$1" -r "$REGION" -a FORMAT/AD -f "$REF" -b "$BAMS" 2>/dev/null \
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
        printf "%s\t%s>%s\tQUAL=%s\t%-9s\t%s\n", $1,$2,$3,$4,verdict,cs
      }' > "$2"
}

dump 1000  "$TMP/d1000"
dump 50000 "$TMP/d50000"

echo
echo "==================================================================="
echo "The 7 flipped sites — REF>ALT, QUAL, verdict, total reads per allele"
echo "(read counts are summed over all 72 samples; 'kept' = passes the"
echo " biallelic + QUAL>59 filter that builds RefAlt.txt)"
echo "==================================================================="
for p in $POSITIONS; do
  a=$(awk -v p="$p" '$1==p{sub(/^[0-9]+\t/,""); print}' "$TMP/d1000")
  b=$(awk -v p="$p" '$1==p{sub(/^[0-9]+\t/,""); print}' "$TMP/d50000")
  echo "chrX:$p"
  printf "   d=1000   %s\n" "${a:-not called here (no variant)}"
  printf "   d=50000  %s\n" "${b:-not called here (no variant)}"
  echo
done
