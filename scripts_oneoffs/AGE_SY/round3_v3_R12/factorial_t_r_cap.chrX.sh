#!/bin/bash
#SBATCH --job-name=factorial_tr
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=10:00:00
#SBATCH -o process/factorial_chrX/task_%a.out
###############################################################################
# factorial_t_r_cap.chrX.sh — 2x2 (flag x cap) over one chrX window, PARALLEL.
#
# The four conditions (-t/-r x -d 1000/50000) run as a 4-task ARRAY (all at
# once), because the -t arms stream the whole BAM and are slow; running them
# sequentially wastes wall clock. A dependent merge job then prints the 2x2.
#
# Submit once from repo root on the cluster (NOT via sbatch — just bash it):
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/factorial_t_r_cap.chrX.sh
# It self-submits the array + the merge. Read the 2x2 in process/factorial_chrX/merge.out
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true

REGION=chrX:15900000-16400000
REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
OUT=process/factorial_chrX
EXAMPLES="16018400 16026914 16028272 16029051 16030427 16030958 16031180 16031709 16034643 16045058"
# condition index -> "flag depth"
COND1="-t 1000"; COND2="-r 1000"; COND3="-t 50000"; COND4="-r 50000"

# ======================= mode: orchestrator (login node) =====================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" && "${1:-}" != "--merge" ]]; then
  SELF=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
  mkdir -p "$OUT"
  JID=$(sbatch --parsable --array=1-4 "$SELF")
  echo "4 conditions submitted as array (parallel): $JID"
  MID=$(sbatch --parsable --dependency=afterok:${JID} \
        -A tdlong_lab -p standard --time=20:00 -o "$OUT/merge.out" \
        --wrap="bash '$SELF' --merge")
  echo "merge (2x2 table) submitted:  $MID"
  echo
  echo "watch:   squeue -u \$USER"
  echo "result:  $OUT/merge.out   (per-arm timings in $OUT/time.*.txt)"
  exit 0
fi

# ============================ mode: merge (verdict) ==========================
if [[ "${1:-}" == "--merge" ]]; then
  echo "==================================================================="
  echo "FACTORIAL RESULT  —  region $REGION"
  echo "Is the old-vs-new caller difference the FETCH FLAG (-t vs -r)"
  echo "or the DEPTH CAP (-d 1000 vs -d 50000)?  Answer below."
  echo "==================================================================="
  echo
  echo "Wall time per condition (this is why -t is the slow one):"
  cat "$OUT"/time.*.txt 2>/dev/null | sed 's/^/  /'
  echo
  awk -v A="$OUT/cond.1.txt" -v B="$OUT/cond.2.txt" -v C="$OUT/cond.3.txt" -v D="$OUT/cond.4.txt" '
    function load(f,tag){ while((getline l < f)>0){ split(l,x,"\t"); pos=x[1]
        na[tag,pos]=x[2]; pass[tag,pos]=x[4]; seen[pos]=1 } close(f) }
    function P(tag,pos){ return ((tag,pos) in pass) ? pass[tag,pos]+0 : 0 }
    function YN(tag,pos){ return ((tag,pos) in pass) ? (pass[tag,pos]+0 ? "SNP":"dropped") : "no-call" }
    BEGIN{
      load(A,"t1"); load(B,"r1"); load(C,"t5"); load(D,"r5")
      for(pos in seen){
        if(P("t1",pos)!=P("r1",pos)) f1++
        if(P("t5",pos)!=P("r5",pos)) f5++
        if(P("t1",pos)!=P("t5",pos) || P("r1",pos)!=P("r5",pos)){ capn++; flip[pos]=1 }
      }
      print "FETCH FLAG  (-t  vs  -r):"
      printf "  positions where -t and -r give a different call at d=1000 : %d\n", f1+0
      printf "  positions where -t and -r give a different call at d=50000: %d\n", f5+0
      print (f1+f5==0) ? "  => the fetch flag changes ZERO calls. It only affects SPEED." \
                       : "  => the fetch flag DOES change calls."
      print ""
      print "DEPTH CAP  (d=1000  vs  d=50000):"
      printf "  positions whose SNP call changes when the cap is lifted: %d\n", capn+0
      if(capn==0) print "  => the cap changes ZERO calls here."
      else {
        print "  => the CAP is responsible. The positions it flips:"
        for(pos in flip)
          printf "     %-10s  d1000: -t=%-7s -r=%-7s | d50000: -t=%-7s -r=%-7s | #alleles capped=%s uncapped=%s\n",
            pos, YN("t1",pos),YN("r1",pos),YN("t5",pos),YN("r5",pos), na["t1",pos], na["t5",pos]
      }
    }
  '
  echo
  echo "key: SNP = passes biallelic + QUAL>59 | dropped = called but filtered out | no-call = not variant"
  exit 0
fi

# ========================= mode: array task (one arm) ========================
i="$SLURM_ARRAY_TASK_ID"
eval "cond=\$COND$i"
read -r flag depth <<< "$cond"
echo "condition $i: mpileup $flag $REGION -d $depth" >&2
t0=$SECONDS
bcftools mpileup -I -d "$depth" "$flag" "$REGION" -a FORMAT/AD -f "$REF" -b "$BAMS" 2>/dev/null \
  | bcftools call -mv 2>/dev/null \
  | bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\n' \
  | awk '{ nalt=($3=="."?0:split($3,a,","))
           issnp=(length($2)==1); for(k=1;k<=nalt;k++) if(length(a[k])!=1) issnp=0
           pass=(nalt==1 && issnp==1 && $4>59)?1:0
           printf "%s\t%d\t%s\t%d\n",$1,nalt,$4,pass }' > "$OUT/cond.$i.txt"
printf "cond %d (%s @ %s): %d s\n" "$i" "$flag" "$depth" "$((SECONDS - t0))" > "$OUT/time.$i.txt"
