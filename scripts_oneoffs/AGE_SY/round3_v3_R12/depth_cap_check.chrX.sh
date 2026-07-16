#!/bin/bash
#SBATCH --job-name=depth_cap_check
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:30:00
#SBATCH -o process/prove_flips/depthcheck_%j.out
###############################################################################
# depth_cap_check.chrX.sh
#
# Direct test of the hypothesis: a flip happens because some sample's TRUE read
# depth exceeds the -d 1000 cap, so mpileup subsamples it, and which reads
# survive moves QUAL. FORMAT/DP is HIGH-QUALITY depth and can look unchanged
# even when a sample is capped, so it can't answer this — we need the true
# (uncapped) read count per sample.
#
# For each of the 7 flip sites (small ±10 kb window, fast):
#   - QUAL from the exact caller at d=1000 and d=50000
#   - TRUE per-sample read depth (samtools, -d 100000 = uncapped)
#   - how many samples are AT/OVER 1000 (i.e. actually capped at d=1000)
#   - for each capped sample: its true depth and its caller DP:AD at both depths
#
# Read literally:
#   samples>=1000 == 0  -> the cap does NOT bind here; the QUAL change is NOT
#                          subsampling -> hypothesis dead for that site.
#   samples>=1000  > 0  -> the cap binds; subsampling is live -> hypothesis holds.
#
# Submit from repo root on the cluster:
#   mkdir -p process/prove_flips
#   sbatch scripts_oneoffs/AGE_SY/round3_v3_R12/depth_cap_check.chrX.sh
###############################################################################
set -uo pipefail
# only bcftools/1.21 in the main shell; samtools/1.10 is isolated in a subshell
# (loading it here would drag in htslib 1.10 and break bcftools).
module load bcftools/1.21 2>/dev/null || true

REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
OUT=process/prove_flips
mkdir -p "$OUT"
POS="16025382 16030427 16030958 16031180 16031851 16045058 16363526"
REGIONS="chrX:16015382-16055058,chrX:16353526-16373526"   # ±10 kb around the 7 sites
CAP=1000

for f in "$REF" "$BAMS"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

# ---- caller QUAL at each depth on the small window ----
for d in 1000 50000; do
  vcf="$OUT/dc.d$d.vcf.gz"
  echo "[$(date +%H:%M:%S)] call -d $d ..." >&2
  bcftools mpileup -I -d "$d" -r "$REGIONS" -a FORMAT/AD,FORMAT/DP -f "$REF" -b "$BAMS" \
    | bcftools call -mv -Oz -o "$vcf"
  bcftools index -f "$vcf"
  n=$(bcftools view -H "$vcf" 2>/dev/null | wc -l | tr -d ' ')
  echo "  -> $n records" >&2
  [[ "$n" -eq 0 ]] && { echo "ERROR: 0 records at d=$d (tool failure). Aborting." >&2; exit 1; }
done
bcftools query -l "$OUT/dc.d1000.vcf.gz" > "$OUT/samples.txt"

qual() { bcftools query -r "chrX:$1-$1" -f '%ALT\t%QUAL\n' "$OUT/dc.d$2.vcf.gz" 2>/dev/null | head -1; }
dpad() { bcftools query -r "chrX:$1-$1" -f '[%SAMPLE\t%DP\t%AD\n]' "$OUT/dc.d$2.vcf.gz" 2>/dev/null \
         | awk -v s="$3" '$1==s{print $2":"$3; exit}'; }

echo "==================================================================="
echo "DIRECT CAP TEST — is any sample's TRUE depth over the $CAP-read cap?"
echo "==================================================================="
for p in $POS; do
  a1=$(qual "$p" 1000);  alt1=${a1%%$'\t'*}; q1=${a1##*$'\t'}
  a2=$(qual "$p" 50000); alt2=${a2%%$'\t'*}; q2=${a2##*$'\t'}
  [[ -z "$a1" ]] && { alt1="-"; q1="no-call"; }
  [[ -z "$a2" ]] && { alt2="-"; q2="no-call"; }
  echo "###################################################################"
  printf "chrX:%s   ALT d1000=%s d50000=%s   QUAL d1000=%s  d50000=%s\n" "$p" "$alt1" "$alt2" "$q1" "$q2"

  # true per-sample read depth, uncapped (samtools isolated in a subshell so its
  # old htslib can't break bcftools in the parent). Capture once, awk it twice.
  raw="$OUT/.raw.$p"
  ( module purge 2>/dev/null; module load samtools/1.10 2>/dev/null
    samtools mpileup -Q0 -B -d 100000 -f "$REF" -r "chrX:$p-$p" -b "$BAMS" 2>/dev/null ) > "$raw"

  read -r maxd over < <(awk -v names="$OUT/samples.txt" -v cap="$CAP" '
      BEGIN{ ni=0; while((getline nm<names)>0) ni++ }
      { max=0; over=0
        for(s=1;s<=ni;s++){ d=$(3*s+1)+0; if(d>max)max=d; if(d>=cap)over++ }
        print max, over }' "$raw")
  : "${maxd:=0}" "${over:=0}"
  printf "  true depth: max=%s   samples >= %s (capped at d=1000): %s\n" "$maxd" "$CAP" "$over"

  if [[ "$over" -gt 0 ]]; then
    awk -v names="$OUT/samples.txt" -v cap="$CAP" '
      BEGIN{ ni=0; while((getline nm<names)>0) name[++ni]=nm }
      { for(s=1;s<=ni;s++){ d=$(3*s+1)+0; if(d>=cap) print name[s]"\t"d } }' "$raw" > "$OUT/.capped.$p"
    printf "    %-20s %-8s %-16s %-16s\n" "sample" "trueDP" "d1000 DP:AD" "d50000 DP:AD"
    while IFS=$'\t' read -r s td; do
      [[ -z "$s" ]] && continue
      printf "    %-20s %-8s %-16s %-16s\n" "$s" "$td" "$(dpad "$p" 1000 "$s")" "$(dpad "$p" 50000 "$s")"
    done < "$OUT/.capped.$p"
    rm -f "$OUT/.capped.$p"
  else
    echo "    -> NO sample hits the $CAP cap here. The cap does not bind; a QUAL change is NOT subsampling."
  fi
  rm -f "$raw"
  echo
done
echo "done."
