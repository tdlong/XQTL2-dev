#!/bin/bash
#SBATCH --job-name=loo_genome
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=4:00:00
#SBATCH -o logs/loo_genome_%j.out
#SBATCH -e logs/loo_genome_%j.out

# Submit from the XQTL2-dev root:  sbatch temp_aging/analysis/loo_genome.sh
#
# 24 GB total. Each meansBySample is ~257 MB on disk and is read one at a time,
# reduced, then dropped -- the four are never in memory together.
#
# Returns its own results through git, as reproduce.sh does, so a finished run
# needs no staging, committing or pushing by hand. Adds ONLY the two output
# files and this job's log -- never `git add -A`, which would sweep up whatever
# else happens to be dirty in the tree. NO_GIT=1 to skip.

module load R/4.2.2

fail=0
Rscript temp_aging/analysis/loo_genome.R || fail=1

if [ "${NO_GIT:-0}" != "1" ]; then
  echo
  echo "── returning results through git ─────────────────────────────"
  git add temp_aging/numbers/loo_positions.txt.gz \
          temp_aging/numbers/loo_tiles.txt.gz 2>/dev/null
  lg="logs/loo_genome_${SLURM_JOB_ID:-manual}.out"
  [ -f "$lg" ] && git add -f "$lg" 2>/dev/null
  if git diff --cached --quiet; then
    echo "   nothing changed"
  else
    git commit -q -m "loo: genome-wide leave-one-cage-out output (job ${SLURM_JOB_ID:-manual})" \
      && echo "   committed" || echo "   commit failed" >&2
    if git pull --rebase -q && git push -q; then
      echo "   pushed"
    else
      echo "   PUSH FAILED -- committed locally, push by hand" >&2
      fail=1
    fi
  fi
fi
exit "$fail"
