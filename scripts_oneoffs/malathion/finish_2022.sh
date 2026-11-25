#!/bin/bash
#SBATCH --job-name=mal_finish
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1:00:00
#SBATCH -o logs/malathion/finish_%j.out
#SBATCH -e logs/malathion/finish_%j.out

# finish_2022.sh — last job of rerun_2022.sh. Runs the paper's ANOVA on the
# all-pools haplotypes, writes the peak tables and figures into logs/malathion/,
# and returns them through git (the only route off the cluster). Submitted by
# rerun_2022.sh with a dependency on the five scans; can also be run by hand:
#   sbatch scripts_oneoffs/malathion/finish_2022.sh <dir> <prefix> <R1 scan,R2 scan,...>
#   e.g.  process/malathion_rep1 malathion_rep1 malathion_rep1_R1,malathion_rep1_R2
#
# Same git discipline as temp_aging/analysis/drift_tiles.sh: nothing is written
# to this log after it is staged, and what git prints goes to an untracked side
# file so the committed log is complete. NO_GIT=1 skips the commit.

module load R/4.2.2
DIR=${1:?project dir}
PREFIX=${2:?scan prefix}
ALSO=${3:-}
OUT=logs/malathion
mkdir -p "$OUT"

fail=0
Rscript scripts_oneoffs/malathion/anova_scan.R \
    --dir  "$DIR" \
    --scan "${PREFIX}_all" \
    --also "$ALSO" \
    --out  "$OUT" || fail=1

SCANS="${PREFIX}_all${ALSO:+,$ALSO}"
for s in ${SCANS//,/ }; do
  f="$DIR/Scans/$s/$s.Manhattan.png"
  if [ -f "$f" ]; then cp "$f" "$OUT/"; else echo "missing $f"; fail=1; fi
done
echo "smoothing R2 per scan:"
for s in ${SCANS//,/ }; do
  printf '  %-22s %s\n' "$s" "$(cat "$DIR/Scans/$s/$s.smooth_r2.txt" 2>/dev/null || echo '?')"
done

if [ "${NO_GIT:-0}" != "1" ]; then
  echo
  lg="$OUT/finish_${SLURM_JOB_ID:-manual}.out"
  gl="$OUT/finish_${SLURM_JOB_ID:-manual}_git.txt"
  echo "── returning results through git; what git says goes to $gl, not here ──"
  exec >"$gl" 2>&1
  git add "$OUT"/*.png "$OUT"/*.txt "$OUT"/*.txt.gz 2>/dev/null
  [ -f "$lg" ] && git add -f "$lg"
  git reset -q -- "$gl" 2>/dev/null
  git status --short
  if git diff --cached --quiet; then
    echo "nothing staged -- nothing to commit"
  else
    git commit -m "malathion 2022 rerun: paper ANOVA vs pipeline Wald, per-pair scans (job ${SLURM_JOB_ID:-manual})" && git push
  fi
fi
exit "$fail"
