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
# Returns the output and this job's SLURM log through git. The log is committed
# so it can be read on the laptop, acted on, then `git rm`'d. Nothing may be
# written to the log after it is added: loo_genome.sh echoed "committed" and
# "pushed" afterwards, so the tracked copy never matched the file and every
# pull on HPC3 failed on it. All output stops at the git add. NO_GIT=1 to skip.

module load R/4.2.2

fail=0
Rscript temp_aging/analysis/drift_tiles.R || fail=1

if [ "${NO_GIT:-0}" != "1" ]; then
  echo
  lg="logs/drift_tiles_${SLURM_JOB_ID:-manual}.out"
  gl="logs/drift_tiles_${SLURM_JOB_ID:-manual}_git.txt"
  echo "── returning results through git; what git says goes to $gl, not here ──"
  # Everything git prints goes to a side file that is never committed, so the
  # committed log is complete and the side file (untracked) cannot block a pull.
  exec >"$gl" 2>&1
  git add temp_aging/numbers/drift_tiles.txt.gz
  [ -f "$lg" ] && git add -f "$lg"
  git status --short
  if git diff --cached --quiet; then
    echo "nothing staged -- nothing to commit"
  else
    git commit -m "drift_tiles: per-replicate pool frequencies on a 0.25 cM grid (job ${SLURM_JOB_ID:-manual})" && git push
  fi
fi
exit "$fail"
