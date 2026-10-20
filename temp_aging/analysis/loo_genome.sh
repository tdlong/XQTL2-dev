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

# TILE_MB=1 tiles in physical distance and writes *_mb files; default is 1 cM.
#   sbatch temp_aging/analysis/loo_genome.sh
#   TILE_MB=1 sbatch --export=ALL,TILE_MB=1 temp_aging/analysis/loo_genome.sh

module load R/4.2.2

fail=0
Rscript temp_aging/analysis/loo_genome.R || fail=1

if [ "${NO_GIT:-0}" != "1" ]; then
  echo
  echo "── returning results through git ─────────────────────────────"
  suf=""; [ "${TILE_MB:-0}" = "1" ] && suf="_mb"
  git add temp_aging/numbers/loo_positions${suf}.txt.gz \
          temp_aging/numbers/loo_tiles${suf}.txt.gz \
          temp_aging/numbers/loo_u${suf}.txt.gz 2>/dev/null
  lg="logs/loo_genome_${SLURM_JOB_ID:-manual}.out"
  [ -f "$lg" ] && git add -f "$lg" 2>/dev/null
  if git diff --cached --quiet; then
    echo "   nothing changed"
  else
    git commit -q -m "loo: genome-wide leave-one-cage-out output (job ${SLURM_JOB_ID:-manual})" \
      && echo "   committed" || echo "   commit failed" >&2
    # No `git pull --rebase` here. It fails on ANY unstaged change, and a job
    # that writes a log always leaves one, so the pull broke every push. Push
    # directly; if origin has genuinely moved, the push says so and Tony
    # resolves it at the login node where he can see the tree.
    if git push -q; then
      echo "   pushed"
    else
      echo "   not pushed -- committed locally, run: git push" >&2
    fi
  fi
fi
exit "$fail"
