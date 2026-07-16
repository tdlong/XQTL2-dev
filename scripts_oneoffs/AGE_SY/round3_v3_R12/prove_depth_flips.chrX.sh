#!/bin/bash
#SBATCH --job-name=prove_depth_flips
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00
#SBATCH -o process/prove_flips/prove_%j.out
###############################################################################
# prove_depth_flips.chrX.sh
#
# QUESTION (settled method): with -r fixed, changing only -d (1000 vs 50000)
# keeps some SNPs and drops others. WHY.
#
# This does not speculate. It runs the EXACT production caller, changing ONLY
# -d, and for every SNP whose kept/dropped verdict flips it prints, side by
# side, the caller's own output at both depths:
#
#   the pipeline caller (bam2bcf2REFALT.tiled.sh:54):
#     bcftools mpileup -I -d D -r REGION -a FORMAT/AD,FORMAT/DP -f ref -b bams
#       | bcftools call -mv
#   the pipeline filter (bam2bcf2REFALT.tiled.sh:63):
#     bcftools view -m2 -M2 -v snps -i 'QUAL>59'
#
# For each flip it reports which of the filter's gates moved, purely from the
# numbers the caller emitted:
#   GATE A (3rd allele) : nALT changed (1<->2). A 3rd allele makes the site
#                         non-biallelic => -m2 -M2 drops it.
#   GATE B (QUAL)       : nALT stayed 1, QUAL crossed 59.
#   GATE C (no-call)    : variant called at one depth, absent at the other.
# Then it prints, for the samples whose depth actually changed under the cap,
# their per-allele read counts (AD) from the caller at both depths, and raw
# per-base pileup counts (samtools, independent of the caller) so you can see
# the exact reads that appeared/vanished.
#
# Submit from repo root on the cluster:
#   mkdir -p process/prove_flips
#   sbatch scripts_oneoffs/AGE_SY/round3_v3_R12/prove_depth_flips.chrX.sh
###############################################################################
set -uo pipefail
# bcftools/1.21 needs htslib>=1.16. Do NOT module-load samtools/1.10 alongside it:
# that drags in htslib/1.10.2 and breaks bcftools (undefined symbol
# bcf_has_variant_types, version HTSLIB_1.16). samtools is loaded ONLY inside an
# isolated subshell for the raw-pileup step, so it can't corrupt bcftools here.
module load bcftools/1.21 2>/dev/null || true

REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
REGION=chrX:15900000-16400000      # contiguous -r window (same one the factorial used)
D1=1000
D2=50000
OUT=process/prove_flips
mkdir -p "$OUT"

for f in "$REF" "$BAMS"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

# ---------------------------------------------------------------------------
# 1) EXACT pipeline call at each depth. Keep the VCFs; every claim below is
#    mined from them.
# ---------------------------------------------------------------------------
for d in "$D1" "$D2"; do
  vcf="$OUT/call.d$d.vcf.gz"
  echo "[$(date +%H:%M:%S)] mpileup|call  -r $REGION  -d $d ..." >&2
  # stderr left visible on purpose: a tool/module failure must be seen, not hidden.
  bcftools mpileup -I -d "$d" -r "$REGION" -a FORMAT/AD,FORMAT/DP -f "$REF" -b "$BAMS" \
    | bcftools call -mv -Oz -o "$vcf"
  bcftools index -f "$vcf"
  n=$(bcftools view -H "$vcf" 2>/dev/null | wc -l | tr -d ' ')
  echo "  -> $n variant records at d=$d" >&2
  if [[ "$n" -eq 0 ]]; then
    echo "ERROR: caller emitted 0 records at d=$d. That is a tool/module failure," >&2
    echo "       not biology (check the bcftools/htslib errors above). Aborting." >&2
    exit 1
  fi
done

bcftools query -l "$OUT/call.d$D1.vcf.gz" > "$OUT/samples.txt"

# ---------------------------------------------------------------------------
# 2) Positions that PASS the exact pipeline filter at each depth; flips =
#    passing at exactly one depth.
# ---------------------------------------------------------------------------
for d in "$D1" "$D2"; do
  bcftools view -m2 -M2 -v snps -i 'QUAL>59' "$OUT/call.d$d.vcf.gz" 2>/dev/null \
    | bcftools query -f '%POS\n' | sort -un > "$OUT/pass.d$d.txt"
done
comm -3 "$OUT/pass.d$D1.txt" "$OUT/pass.d$D2.txt" | tr -d '\t' | sort -un > "$OUT/flips.pos"

NFLIP=$(wc -l < "$OUT/flips.pos")
echo "==================================================================="
echo "REGION $REGION   caller: bcftools mpileup -I | call -mv"
echo "filter: -m2 -M2 -v snps -i 'QUAL>59'    (only -d differs: $D1 vs $D2)"
echo "kept @ d=$D1: $(wc -l < "$OUT/pass.d$D1.txt")   kept @ d=$D2: $(wc -l < "$OUT/pass.d$D2.txt")"
echo "FLIPS (kept at exactly one depth): $NFLIP"
echo "==================================================================="
echo

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
# site record at one pos in one vcf -> "REF<TAB>ALT<TAB>nALT<TAB>QUAL" (SNP rec), or empty
rec() { # $1=pos $2=vcf
  bcftools query -r "chrX:$1-$1" -f '%REF\t%ALT\t%QUAL\n' "$2" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"} { if(length($1)!=1) next
             na=($2=="."?0:split($2,a,",")); snp=1; for(i=1;i<=na;i++) if(length(a[i])!=1) snp=0
             if(snp){ print $1,$2,na,$3; exit } }'
}
verdict() { # $1=nALT $2=QUAL  -> KEPT / DROPPED / (caller decides biallelic+snp already ensured by rec)
  awk -v n="$1" -v q="$2" 'BEGIN{ print (n==1 && q+0>59) ? "KEPT" : "DROPPED" }'
}

# ---------------------------------------------------------------------------
# 3) per flip: the proof
# ---------------------------------------------------------------------------
while read -r pos; do
  [[ -z "$pos" ]] && continue
  r1=$(rec "$pos" "$OUT/call.d$D1.vcf.gz"); r2=$(rec "$pos" "$OUT/call.d$D2.vcf.gz")
  IFS=$'\t' read -r ref1 alt1 na1 q1 <<< "$r1"
  IFS=$'\t' read -r ref2 alt2 na2 q2 <<< "$r2"
  [[ -z "${na1:-}" ]] && { na1=0; alt1="-"; q1="-"; ref1="${ref2:-?}"; }
  [[ -z "${na2:-}" ]] && { na2=0; alt2="-"; q2="-"; ref2="${ref1:-?}"; }
  v1="no-call"; [[ "$na1" -ge 1 ]] && v1=$(verdict "$na1" "$q1")
  v2="no-call"; [[ "$na2" -ge 1 ]] && v2=$(verdict "$na2" "$q2")

  # which gate (mechanical, from the two records)
  if [[ "$na1" -eq 0 || "$na2" -eq 0 ]]; then gate="GATE C: no-call at one depth (variant not called)"
  elif [[ "$na1" -ne "$na2" ]]; then gate="GATE A: 3rd allele — nALT $na1 vs $na2 => -m2 -M2 drops the multiallelic side"
  elif [[ "$v1" != "$v2" ]]; then gate="GATE B: QUAL crossed 59 ($q1 vs $q2), nALT stayed 1"
  else gate="(verdict same at both depths — not a flip; check filter)"
  fi

  echo "###################################################################"
  printf "chrX:%s   REF=%s\n" "$pos" "${ref2:-$ref1}"
  printf "  d=%-6s ALT=%-10s nALT=%s QUAL=%-9s -> %s\n" "$D1" "$alt1" "$na1" "$q1" "$v1"
  printf "  d=%-6s ALT=%-10s nALT=%s QUAL=%-9s -> %s\n" "$D2" "$alt2" "$na2" "$q2" "$v2"
  printf "  >>> %s\n" "$gate"

  # per-sample AD from the caller, at both depths, for samples whose DP changed
  # (those are the samples the cap actually touched — the cause).
  bcftools query -r "chrX:$pos-$pos" -f '[%SAMPLE\t%DP\t%AD\n]' "$OUT/call.d$D1.vcf.gz" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' > "$OUT/.ad1"
  bcftools query -r "chrX:$pos-$pos" -f '[%SAMPLE\t%DP\t%AD\n]' "$OUT/call.d$D2.vcf.gz" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' > "$OUT/.ad2"
  echo "  samples the cap changed (DP differs) — AD is per-allele in the header's allele order:"
  printf "    %-20s %-22s %-22s\n" "sample" "d$D1 DP:AD" "d$D2 DP:AD"
  join -t $'\t' -1 1 -2 1 <(sort -k1,1 "$OUT/.ad1") <(sort -k1,1 "$OUT/.ad2") \
    | awk -F'\t' -v OFS='' '$2!=$4 {printf "    %-20s %-22s %-22s\n", $1, $2":"$3, $4":"$5}'

  # raw base counts, independent of the caller (samtools), for those same samples.
  # samtools/1.10 is loaded in an ISOLATED subshell (module purge first) so its
  # old htslib can never leak back and break bcftools in the parent shell.
  capped=$(join -t $'\t' -1 1 -2 1 <(sort -k1,1 "$OUT/.ad1") <(sort -k1,1 "$OUT/.ad2") \
           | awk -F'\t' '$2!=$4 {print $1}')
  if [[ -n "$capped" ]]; then
    echo "  raw pileup base counts (samtools, min-BQ0, no BAQ — independent of the caller):"
    for d in "$D1" "$D2"; do
      ( module purge 2>/dev/null; module load samtools/1.10 2>/dev/null
        samtools mpileup -Q0 -B -d "$d" -f "$REF" -r "chrX:$pos-$pos" -b "$BAMS" 2>/dev/null ) \
          | awk -v names="$OUT/samples.txt" -v cap="$capped" -v dd="$d" '
              function strip(s,   out,i,c,j,num,l){
                gsub(/\^./,"",s); gsub(/\$/,"",s)              # read-start (^ + mapq) and read-end
                out=""; i=1
                while(i<=length(s)){ c=substr(s,i,1)
                  if(c=="+"||c=="-"){ j=i+1; num=""
                      while(j<=length(s) && substr(s,j,1)~/[0-9]/){num=num substr(s,j,1);j++}
                      l=num+0; i=j+l }                          # skip indel: +/- , digits, l bases
                  else { out=out c; i++ } }
                return out }
              function cnt(s,ref,   t,A,C,G,T,r,i,ch){
                t=strip(s); r=gsub(/[.,]/,"",t)
                A=gsub(/[Aa]/,"",t); C=gsub(/[Cc]/,"",t); G=gsub(/[Gg]/,"",t); T=gsub(/[Tt]/,"",t)
                return sprintf("REF(%s):%d A:%d C:%d G:%d T:%d", ref, r, A, C, G, T) }
              BEGIN{ ni=0; while((getline nm<names)>0){ nm2[++ni]=nm }
                     split(cap,cc,"\n"); for(k in cc) want[cc[k]]=1 }
              { ref=$3
                for(s=1;s<=ni;s++){ dp=$(3*s+1); bs=$(3*s+2)
                  if((nm2[s] in want)) printf "    d%-6s %-20s DP=%-6s %s\n", dd, nm2[s], dp, cnt(bs,ref) } }'
    done
  fi
  echo
done < "$OUT/flips.pos"

rm -f "$OUT/.ad1" "$OUT/.ad2"
echo "done. raw VCFs kept at $OUT/call.d$D1.vcf.gz and $OUT/call.d$D2.vcf.gz"
