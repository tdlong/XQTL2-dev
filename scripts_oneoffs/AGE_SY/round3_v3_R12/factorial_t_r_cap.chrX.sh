#!/bin/bash
#SBATCH --job-name=factorial_tr
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=8:00:00
#SBATCH -o factorial_t_r_cap.chrX.%j.out
###############################################################################
# factorial_t_r_cap.chrX.sh — 2x2 experiment over one 500 kb chrX window.
#
# Run bcftools mpileup|call FOUR ways and compare the SNP calls:
#     fetch flag:  -t  vs  -r
#     max depth:   -d 1000  vs  -d 50000
#
# The window chrX:15,900,000-16,400,000 contains ~10 sites where the whole-genome
# run had the two callers disagree. For each variant position we record, per
# condition: number of ALT alleles, QUAL, and whether it PASSES the RefAlt filter
# (biallelic SNP, QUAL>59). Then we show every position where the four conditions
# do not all agree.
#
# Expected reads:
#   -t@50000 == -r@50000 (uncapped: flag can't matter)  -> flag is innocent
#   calls change between @1000 and @50000 at the known sites -> the cap is the cause
#
# Submit from repo root on the cluster:
#   sbatch scripts_oneoffs/AGE_SY/round3_v3_R12/factorial_t_r_cap.chrX.sh
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true

CHR=chrX
REGION=chrX:15900000-16400000
REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
# known disagreeing sites inside the window (from the whole-genome run)
EXAMPLES="16018400 16026914 16028272 16029051 16030427 16030958 16031180 16031709 16034643 16045058"

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# one condition -> file of: POS  nALT  QUAL  PASS(1/0)
run () {  # $1 = flag (-t|-r)   $2 = depth   $3 = outfile
  echo "  running mpileup $1 $REGION -d $2 ..." >&2
  bcftools mpileup -I -d "$2" "$1" "$REGION" -a FORMAT/AD -f "$REF" -b "$BAMS" 2>/dev/null \
    | bcftools call -mv 2>/dev/null \
    | bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\n' \
    | awk '{ nalt=($3=="."?0:split($3,a,","))
             issnp=(length($2)==1); for(i=1;i<=nalt;i++) if(length(a[i])!=1) issnp=0
             pass=(nalt==1 && issnp==1 && $4>59)?1:0
             printf "%s\t%d\t%s\t%d\n",$1,nalt,$4,pass }' > "$3"
}

run -t 1000  "$TMP/t.1k"
run -r 1000  "$TMP/r.1k"
run -t 50000 "$TMP/t.50k"
run -r 50000 "$TMP/r.50k"

echo
echo "=================================================================="
echo "2x2 over $REGION   (each cell: nALT/QUAL/PASS)"
echo "=================================================================="
printf "%-10s  %-14s %-14s %-14s %-14s  %s\n" "POS" "-t @1000" "-r @1000" "-t @50000" "-r @50000" "note"
awk '
  function cell(tag,pos){ if(!((tag,pos) in na)) return "----"; return na[tag,pos]"/"q[tag,pos]"/"(p[tag,pos]?"YES":"no") }
  function load(f,tag){ while((getline l < f)>0){ split(l,x,"\t"); pos=x[1]
      na[tag,pos]=x[2]; q[tag,pos]=x[3]; p[tag,pos]=x[4]; seen[pos]=1 } close(f) }
  BEGIN{
    load(A,"t1"); load(B,"r1"); load(C,"t5"); load(D,"r5")
    ne=split(EX,ev," "); for(i=1;i<=ne;i++) isex[ev[i]]=1
    for(pos in seen){
      differ = !(p["t1",pos]==p["r1",pos] && p["r1",pos]==p["t5",pos] && p["t5",pos]==p["r5",pos])
      if(differ || (pos in isex)){
        note=(pos in isex)?"<-- known disagreement":""
        if(differ) note=note "  [conditions differ]"
        printf "%-10s  %-14s %-14s %-14s %-14s  %s\n", pos, cell("t1",pos),cell("r1",pos),cell("t5",pos),cell("r5",pos), note
      }
    }
  }
' A="$TMP/t.1k" B="$TMP/r.1k" C="$TMP/t.50k" D="$TMP/r.50k" EX="$EXAMPLES" | sort -k1,1n

echo
echo "PASS counts per condition (SNPs passing biallelic+QUAL>59 in the window):"
for f in "$TMP/t.1k:-t@1000" "$TMP/r.1k:-r@1000" "$TMP/t.50k:-t@50000" "$TMP/r.50k:-r@50000"; do
  printf "  %-12s %d\n" "${f#*:}" "$(awk '$4==1' "${f%%:*}" | wc -l)"
done
