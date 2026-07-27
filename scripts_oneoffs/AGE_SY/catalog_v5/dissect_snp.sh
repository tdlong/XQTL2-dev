#!/bin/bash
#SBATCH --job-name=dissect_snp
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:20:00
#SBATCH -o logs/AGE_SY/dissect_%j.out
###############################################################################
# dissect_snp.sh — read-by-read forensic of ONE (SNP,sample) where v3 and v4
# disagree. No summaries: dump every read's base / base-quality / MAPQ at the
# position, then re-count under each recipe's flags toggled one at a time, so the
# exact filter that drops reads is visible. Default site: chrX:13605110 Con_R5_F
# (v3 4/76 = 0.05 ALT vs v4 146/218 = 0.67 ALT — v3 is missing ~142 reads).
#
#   CHR=chrX POS=13605110 BAM=data/bam/AGE_SY/Con_R5_F.bam \
#     sbatch scripts_oneoffs/AGE_SY/catalog_v5/dissect_snp.sh
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true
module load samtools/1.21 2>/dev/null || true
REF=pipeline/ref/dm6.fa
CHR=${CHR:-chrX}
POS=${POS:-13605110}
BAM=${BAM:-data/bam/AGE_SY/Con_R5_F.bam}
CAT=${CAT:-process/AGE_SY_v5/Catalog/catalog.tsv.gz}
SAMP=$(basename "$BAM" .bam)
for f in "$REF" "$BAM"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

echo "########## DISSECT $CHR:$POS   sample $SAMP ##########"
RB=$(samtools faidx "$REF" "$CHR:$POS-$POS" | tail -1 | tr a-z A-Z)
echo "REF base (dm6)      : $RB"
echo -n "catalog REF,ALT     : "; tabix "$CAT" "$CHR:$POS-$POS" 2>/dev/null | awk '{print $3}' || echo "(catalog not tabixed)"
echo

echo "=== nearest indel evidence within +/-50 bp (positions where reads show +/-N) ==="
samtools mpileup -B -q0 -Q0 -f "$REF" -r "$CHR:$((POS-50))-$((POS+50))" "$BAM" 2>/dev/null \
  | awk -v p=$POS 'match($5,/[+-][0-9]/){d=$2-p; print "  indel@"$2"  (offset "d")  pileup="substr($5,1,40)}'
echo "  (none printed => no indel within 50 bp)"
echo

echo "=== recipe re-count at the single site — FORMAT/AD (ref,alt reads each recipe COUNTS) ==="
adline() {  # $1=label  $2..=mpileup flags
  local lab="$1"; shift
  local o
  o=$(bcftools mpileup "$@" -r "$CHR:$POS-$POS" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
        | bcftools query -f '%REF\t%ALT\t[%AD]\n' 2>/dev/null | head -1)
  printf "  %-26s REF/ALT=%s  AD(ref,alt)=%s\n" "$lab" "$(echo "$o"|cut -f1-2|tr '\t' '/')" "$(echo "$o"|cut -f3)"
}
adline "raw   (-B -q0  -Q0 )"  -B -q0  -Q0
adline "v3    (BAQ -q0  -Q13)" -I -d 1000 -q0  -Q13
adline "v3-BAQ off (-B -Q13)"  -I -d 1000 -B -q0  -Q13
adline "isolate -Q20 (BAQ off)" -B -q0  -Q20
adline "isolate -q20 (BAQ off)" -B -q20 -Q13
adline "v4    (-B -q20 -Q20)"  -B -q20 -Q20 --max-depth 2000
echo "  read down the ladder: the row where ALT collapses is the filter responsible."
echo

echo "=== per-read table: BASE  baseQ  MAPQ  (no filters, so every read shown) ==="
echo "    (REF base='$RB'; '.'/',' = ref; a letter = that base; split stats below)"
samtools mpileup -B -q0 -Q0 -f "$REF" -r "$CHR:$POS-$POS" -s --output-QNAME "$BAM" 2>/dev/null \
| awk -v RB="$RB" '
function ord(c,  i){return I[c]}
BEGIN{for(i=0;i<256;i++) I[sprintf("%c",i)]=i}
{
  bases=$5; bq=$6; mq=$7;
  nb=0; delete B;
  # parse pileup base string into per-read bases (handle ^ $ +N -N *)
  L=length(bases);
  for(i=1;i<=L;i++){
    c=substr(bases,i,1);
    if(c=="^"){i++; continue}                 # read start: skip the map-qual char
    if(c=="$"){continue}                       # read end marker
    if(c=="+"||c=="-"){                         # indel: +N or -N then N bases
      j=i+1; num="";
      while(substr(bases,j,1) ~ /[0-9]/){num=num substr(bases,j,1); j++}
      i=j+num-1; continue
    }
    nb++; B[nb]=c;                              # a real base column (incl * deletion)
  }
  # walk bases with aligned baseQ / mapQ (one char each per real base)
  refc=0; altc=0; sumMQalt=0; sumBQalt=0; sumMQref=0; sumBQref=0;
  for(k=1;k<=nb;k++){
    ch=B[k]; q=ord(substr(bq,k,1))-33; m=ord(substr(mq,k,1))-33;
    up=toupper(ch);
    isref=(ch=="."||ch==",");
    isalt=(up ~ /[ACGT]/ && up!=RB);
    tag=(isref?"ref":(isalt?"ALT":"oth"));
    if(printed++<400) printf "  %-3s base=%s baseQ=%2d MAPQ=%2d\n", tag, up, q, m;
    if(isref){refc++; sumMQref+=m; sumBQref+=q}
    else if(isalt){altc++; sumMQalt+=m; sumBQalt+=q}
  }
  print "";
  printf "SPLIT  ref reads=%d  meanMAPQ=%.1f meanBaseQ=%.1f\n", refc, (refc?sumMQref/refc:0), (refc?sumBQref/refc:0);
  printf "SPLIT  ALT reads=%d  meanMAPQ=%.1f meanBaseQ=%.1f\n", altc, (altc?sumMQalt/altc:0), (altc?sumBQalt/altc:0);
  printf "SPLIT  total real bases parsed=%d\n", nb;
}'
echo
echo "=== samtools tview ASCII (text mode, window around site) ==="
COLUMNS=200 samtools tview -d T -p "$CHR:$POS" "$BAM" "$REF" 2>/dev/null | head -50
echo "########## END ##########"
