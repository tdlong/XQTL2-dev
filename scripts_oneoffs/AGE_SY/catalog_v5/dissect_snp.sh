#!/bin/bash
#SBATCH --job-name=dissect_snp
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:20:00
#SBATCH -o logs/AGE_SY/dissect_%j.out
###############################################################################
# dissect_snp.sh — STEP 3. One (SNP,sample) where v3 and v4 disagree and the
# founders say the SNP is CLEAN. Show the reads and find the exact filter that
# zeroes v4's ALT. Default: chrX:6471332 Con_R5_F (v3 56/291 ALT, v4 0/227 ALT;
# dist_indel 81; founders 4 REF-fixed / 4 ALT-fixed).
#
# Prints:
#  A. REF base, catalog REF/ALT, and the ALT base actually seen in the reads.
#  B. FLAG LADDER (unconstrained mpileup AD) to isolate -q / -Q / BAQ:
#       raw(-B -q0 -Q0) -> v3(BAQ -q0 -Q13) -> -B -q0 -Q13 -> -Q20 -> -q20 -> v4flags(-B -q20 -Q20)
#  C. REAL callers: v4 FULL (-B -q20 -Q20 -T cat | call -m -C alleles) and v3 FULL (call -mv).
#     If v4flags(B) already zeroes ALT -> a quality filter did it; if only v4 FULL zeroes it
#     -> the -C alleles catalog constraint did it.
#  D. PER-READ table: base / baseQ / MAPQ, split REF-base vs ALT-base (the "why": are the
#     ALT reads low-MAPQ? -> -q20 removes them = reference mapping bias).
#  E. samtools tview ASCII pileup.
#
#   CHR=chrX POS=6471332 BAM=data/bam/AGE_SY/Con_R5_F.bam \
#     sbatch scripts_oneoffs/AGE_SY/catalog_v5/dissect_snp.sh
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true
module load samtools/1.10 2>/dev/null || true
REF=pipeline/ref/dm6.fa
CHR=${CHR:-chrX}
POS=${POS:-6471332}
BAM=${BAM:-data/bam/AGE_SY/Con_R5_F.bam}
CAT=${CAT:-process/AGE_SY_v5/Catalog/catalog.tsv.gz}
SAMP=$(basename "$BAM" .bam)
for f in "$REF" "$BAM"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done
R="$CHR:$POS-$POS"

echo "########## DISSECT $CHR:$POS   sample $SAMP ##########"
RB=$(samtools faidx "$REF" "$R" | tail -1 | tr a-z A-Z)
echo "A. REF base (dm6)   : $RB"
echo -n "   catalog REF,ALT  : "; tabix "$CAT" "$R" 2>/dev/null | awk '{print $3}' || echo "(not found)"
echo "   observed bases (raw pileup, base:count):"
samtools mpileup -B -q0 -Q0 -f "$REF" -r "$R" "$BAM" 2>/dev/null \
 | awk -v RB="$RB" '{n=split(toupper($5),c,""); for(i=1;i<=n;i++){b=c[i];
     if(b=="."||b==","){r++} else if(b~/[ACGT]/){cnt[b]++}}
     printf "     REF(%s)=%d  ", RB, r; for(k in cnt) printf "%s=%d  ",k,cnt[k]; print ""}'
echo

echo "B. FLAG LADDER — unconstrained mpileup FORMAT/AD (REF ALT AD):"
mp(){ local lab="$1"; shift
  local o; o=$(bcftools mpileup "$@" -r "$R" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
        | bcftools query -f '%REF\t%ALT\t[%AD]\n' 2>/dev/null | head -1)
  printf "   %-26s %s\n" "$lab" "$o"; }
mp "raw    -B -q0  -Q0"   -B -q0  -Q0
mp "v3flag  BAQ -q0 -Q13"     -q0  -Q13
mp "        -B -q0  -Q13" -B -q0  -Q13
mp "isolate -Q20 (BAQ off)" -B -q0  -Q20
mp "isolate -q20 (BAQ off)" -B -q20 -Q13
mp "v4flag  -B -q20 -Q20"  -B -q20 -Q20
echo

echo "C. REAL callers (the actual recipes):"
V4=$(bcftools mpileup -B -q20 -Q20 --max-depth 2000 -r "$R" -T "$CAT" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
     | bcftools call -m -C alleles -T "$CAT" 2>/dev/null | bcftools query -f '%REF\t%ALT\t[%AD]\n' | head -1)
V3=$(bcftools mpileup -I -d 1000 -r "$R" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
     | bcftools call -mv 2>/dev/null | bcftools query -f '%REF\t%ALT\t[%AD]\n' | head -1)
printf "   %-26s %s\n" "v4 FULL (-C alleles)" "${V4:-<no record>}"
printf "   %-26s %s\n" "v3 FULL (call -mv)"   "${V3:-<no variant record>}"
echo "   -> if v4flag(B) already zeroes ALT: a quality filter. If only v4 FULL zeroes it: -C alleles."
echo

echo "D. PER-READ (no filters): base  baseQ  MAPQ   [split REF vs ALT below]"
samtools mpileup -B -q0 -Q0 -f "$REF" -r "$R" -s "$BAM" 2>/dev/null \
| awk -v RB="$RB" '
  BEGIN{for(i=0;i<256;i++) O[sprintf("%c",i)]=i}
  {
    bases=$5; bq=$6; mq=$7; nb=0; L=length(bases);
    for(i=1;i<=L;i++){ c=substr(bases,i,1);
      if(c=="^"){i++; continue} if(c=="$") continue;
      if(c=="+"||c=="-"){ j=i+1; num=""; while(substr(bases,j,1)~/[0-9]/){num=num substr(bases,j,1);j++} i=j+num-1; continue }
      nb++; B[nb]=c }
    for(k=1;k<=nb;k++){ ch=B[k]; q=O[substr(bq,k,1)]-33; m=O[substr(mq,k,1)]-33; up=toupper(ch);
      isref=(ch=="."||ch==","); isalt=(up~/[ACGT]/ && up!=RB);
      if(isref){rc++;rq+=q;rm+=m; rmlo+=(m<20)} else if(isalt){ac++;aq+=q;am+=m; amlo+=(m<20);
        if(showa++<60) printf "   ALT  base=%s baseQ=%2d MAPQ=%2d\n",up,q,m } }
    print "";
    printf "   SPLIT REF-base reads: %d  meanBaseQ=%.1f meanMAPQ=%.1f  (MAPQ<20: %d)\n", rc,(rc?rq/rc:0),(rc?rm/rc:0),rmlo;
    printf "   SPLIT ALT-base reads: %d  meanBaseQ=%.1f meanMAPQ=%.1f  (MAPQ<20: %d)\n", ac,(ac?aq/ac:0),(ac?am/ac:0),amlo;
  }'
echo
echo "E. samtools tview:"
COLUMNS=180 samtools tview -d T -p "$CHR:$POS" "$BAM" "$REF" 2>/dev/null | head -40
echo "########## END ##########"