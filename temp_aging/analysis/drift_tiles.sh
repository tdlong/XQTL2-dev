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
# Returns the output and a log through git. The log is the R script's output,
# written to temp_aging/numbers and closed before the commit -- NOT the SLURM
# .out, which SLURM keeps appending to after the commit and so is dirty forever
# once tracked (what loo_genome.sh did). NO_GIT=1 to skip.

module load R/4.2.2

fail=0
LOG=temp_aging/numbers/drift_tiles.log
Rscript temp_aging/analysis/drift_tiles.R > "$LOG" 2>&1 || fail=1
cat "$LOG"

if [ "${NO_GIT:-0}" != "1" ]; then
  echo
  echo "── returning results through git ─────────────────────────────"
  git add temp_aging/numbers/drift_tiles.txt.gz "$LOG" 2>/dev/null
  if git diff --cached --quiet; then
    echo "   nothing changed"
  else
    git commit -q -m "drift_tiles: per-replicate pool frequencies on a 0.25 cM grid (job ${SLURM_JOB_ID:-manual})" \
      && echo "   committed" || echo "   commit failed" >&2
    if git push -q; then echo "   pushed"; else echo "   not pushed -- committed locally, run: git push" >&2; fi
  fi
fi
exit "$fail"
