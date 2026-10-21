#!/bin/bash
#SBATCH --job-name=drift_tiles
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1:00:00
#SBATCH -o logs/drift_tiles_%j.out
#SBATCH -e logs/drift_tiles_%j.out

# Submit from the XQTL2-dev root:  sbatch temp_aging/analysis/drift_tiles.sh
# Same inputs and memory as loo_genome.sh (which ran in about four minutes).
# Returns the output and its log through git, as loo_genome.sh does. NO_GIT=1 to skip.

module load R/4.2.2

fail=0
Rscript temp_aging/analysis/drift_tiles.R || fail=1

if [ "${NO_GIT:-0}" != "1" ]; then
  echo
  echo "── returning results through git ─────────────────────────────"
  git add temp_aging/numbers/drift_tiles.txt.gz 2>/dev/null
  lg="logs/drift_tiles_${SLURM_JOB_ID:-manual}.out"
  [ -f "$lg" ] && git add -f "$lg" 2>/dev/null
  if git diff --cached --quiet; then
    echo "   nothing changed"
  else
    git commit -q -m "drift_tiles: per-replicate pool frequencies on a 0.25 cM grid (job ${SLURM_JOB_ID:-manual})" \
      && echo "   committed" || echo "   commit failed" >&2
    if git push -q; then echo "   pushed"; else echo "   not pushed -- committed locally, run: git push" >&2; fi
  fi
fi
exit "$fail"
