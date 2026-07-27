#!/bin/bash
#SBATCH --job-name=founder_triage
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00
#SBATCH -o logs/AGE_SY/founder_triage_%j.out
###############################################################################
# founder_triage.sh — STEP 2. Take the big v3-vs-v4 disagreements (>10% ALT-freq
# gap) and run each SNP through the NEW founder catalog to split them into:
#   (a) RED-FLAGGED: the founders show a reason it should never have been called
#       (indel <25bp, min founder depth <10, an intermediate founder, or not
#       segregating) -> these are filter recommendations.
#   (b) CLEAN: no founder red flag -> the remainder we take to the alignments (Step 3).
# Also crosses the split with disagreement direction (v4 drops alt vs v3 drops alt).
#
# INPUTS : process/AGE_SY_v4/Calls/delta_gt2.chrX.tsv.gz (cols: chr pos sample total
#          delta ref_v3 alt_v3 ref_v4 alt_v4), process/AGE_SY_v5/Catalog/catalog.annot.tsv.gz
# OUTPUT : printed report (this log) + logs/AGE_SY/big_disagree_clean.chrX.tsv
#          (clean remainder, worst sample per SNP, sorted worst-first -> Step 3 pool).
#
#   sbatch scripts_oneoffs/AGE_SY/catalog_v5/founder_triage.sh
###############################################################################
set -uo pipefail
ANNOT=process/AGE_SY_v5/Catalog/catalog.annot.tsv.gz
DELTA=process/AGE_SY_v4/Calls/delta_gt2.chrX.tsv.gz
CLEAN=logs/AGE_SY/big_disagree_clean.chrX.tsv
for f in "$ANNOT" "$DELTA"; do [[ -s "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

echo "###### STEP 2: founder triage of chrX big disagreements (|ALT-freq gap|>0.10) ######"
echo "red flags: near_indel(dist<25) | low_founder_depth(min<10) | intermediate_founder | not_segregating"
echo

{ zcat "$ANNOT"; echo "@@@DELTA@@@"; zcat "$DELTA"; } | awk -F'\t' -v CLEAN="$CLEAN" '
  BEGIN{ mode=0 }
  $0=="@@@DELTA@@@" { mode=1; next }
  # ---- annot pass: per-pos founder flags ----
  mode==0 {
    if(FNR==1) next;                                   # annot header
    d=$5+0; mind=1e9; nint=0; nref=0; nalt=0;
    for(i=6;i<=NF;i++){ split($i,a,","); t=a[1]+a[2]; if(t<mind)mind=t;
      if(t==0) continue; f=a[2]/t; if(f<=0.03)nref++; else if(f>=0.97)nalt++; else nint++ }
    D[$2]=d; MIND[$2]=mind; NINT[$2]=nint; NREF[$2]=nref; NALT[$2]=nalt; A[$2]=1;
    next
  }
  # ---- delta pass: keep big pairs, track worst sample per pos ----
  $1=="chr" { next }
  {
    t3=$6+$7; t4=$8+$9; if(t3<=0||t4<=0) next;
    f3=$7/t3; f4=$9/t4; dd=f4-f3; ad=(dd<0?-dd:dd);
    if(ad<=0.10) next;
    p=$2; bigpair++;
    if(!(p in BIG)){ BIG[p]=1; uniq++ }
    if(dd>0) v4hi[p]++; else v3hi[p]++;                # v4>v3 => v3 dropped alt
    if(ad>bestad[p]){ bestad[p]=ad; best[p]=sprintf("%s\tv3:%d/%d\tv4:%d/%d\t%s\t%.3f",
        $3,$7,t3,$9,t4,(dd>0?"v4hi":"v3hi"),ad) }
  }
  END{
    print "chr\tpos\tdist_indel\tmin_fdepth\tnRfix\tnAfix\tnInter\tworst_sample\tv3_alt/tot\tv4_alt/tot\tdir\tgap" > CLEAN;
    for(p in BIG){
      if(!(p in A)){ nia++; nia_v4hi+=(v4hi[p]?v4hi[p]:0); nia_v3hi+=(v3hi[p]?v3hi[p]:0); continue }
      d=D[p]; mind=MIND[p]; nint=NINT[p]; nref=NREF[p]; nalt=NALT[p];
      ni=(d<25); ld=(mind<10); it=(nint>0); ns=(nref==0||nalt==0);
      if(ni)cni++; if(ld)cld++; if(it)cit++; if(ns)cns++;
      any=(ni||ld||it||ns);
      if(any){ cany++; redv4+=(v4hi[p]?v4hi[p]:0); redv3+=(v3hi[p]?v3hi[p]:0) }
      else   { cclean++; clnv4+=(v4hi[p]?v4hi[p]:0); clnv3+=(v3hi[p]?v3hi[p]:0);
               print "chrX\t"p"\t"d"\t"mind"\t"nref"\t"nalt"\t"nint"\t"best[p] > CLEAN }
    }
    printf "BIG disagreements: %d (SNP,sample) pairs over %d unique SNPs\n\n", bigpair, uniq;
    print  "Founder triage of the unique SNPs:";
    printf "  near_indel(<25)        : %7d\n", cni;
    printf "  low_founder_depth(<10) : %7d\n", cld;
    printf "  intermediate_founder   : %7d\n", cit;
    printf "  not_segregating        : %7d\n", cns;
    printf "  ------------------------------------\n";
    printf "  ANY red flag           : %7d (%5.1f%% of unique)\n", cany, 100*cany/uniq;
    printf "  CLEAN (no red flag)    : %7d (%5.1f%% of unique)\n", cclean, 100*cclean/uniq;
    printf "  not in v5 annot        : %7d\n", nia;
    printf "\nDirection x founder status (PAIRS):\n";
    printf "  CLEAN    : v4-drops-alt=%d  v3-drops-alt=%d\n", clnv4, clnv3;
    printf "  RED-FLAG : v4-drops-alt=%d  v3-drops-alt=%d\n", redv4, redv3;
    printf "\nclean remainder list -> %s (sort -k12 -rn for worst-first)\n", CLEAN;
  }
'
echo "###### END STEP 2 ######"