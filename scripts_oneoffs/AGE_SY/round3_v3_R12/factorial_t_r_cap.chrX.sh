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

# ============================ mode: merge (2x2) ==============================
if [[ "${1:-}" == "--merge" ]]; then
  echo "per-arm wall times:"; cat "$OUT"/time.*.txt 2>/dev/null | sed 's/^/  /'
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
  ' A="$OUT/cond.1.txt" B="$OUT/cond.2.txt" C="$OUT/cond.3.txt" D="$OUT/cond.4.txt" EX="$EXAMPLES" | sort -k1,1n
  echo
  echo "PASS counts per condition:"
  for i in 1 2 3 4; do printf "  cond %d : %d\n" "$i" "$(awk '$4==1' "$OUT/cond.$i.txt" | wc -l)"; done
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
