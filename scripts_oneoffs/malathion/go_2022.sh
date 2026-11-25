#!/bin/bash
#SBATCH --job-name=mal_driver
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=3-00:00:00
#SBATCH -o logs/malathion/driver_%j.out
#SBATCH -e logs/malathion/driver_%j.out

# go_2022.sh — the whole 2022 malathion run from one submission:
#   sbatch scripts_oneoffs/malathion/go_2022.sh
# Polls every 5 minutes. Resubmits fetch_sra.sh until all 16 reads are verified,
# builds the A-founder catalog if nothing is building it, and when both are
# there runs the standard pipeline (run_full_pipeline.sh). Then exits. Its log
# is logs/malathion/driver_<jobid>.out.
set -uo pipefail
shopt -s nullglob
[ -d scripts_oneoffs/malathion ] || { echo "ERROR: run from the XQTL2-dev root" >&2; exit 1; }
# Jobs submitted from inside a job inherit this job's SLURM_* environment, which
# would override their own resource requests. Drop it.
unset $(env | sed -n 's/^\(SLURM_[^=]*\)=.*/\1/p')

RAW=data/raw/malathion_2022
CATALOG=process/malathion_2022/Catalog
MAX_FETCH=12
n_fetch=0

# leftovers of an earlier chain that can never start
scancel --name=mal_part2 2>/dev/null
scancel --name=bwa 2>/dev/null

while :; do
  ok=("$RAW"/*.md5ok); n_ok=${#ok[@]}
  if [ "$n_ok" -lt 16 ]; then
    if squeue -u "$USER" -h -n sra_fetch | grep -q .; then
      echo "$(date '+%F %T')  reads $n_ok/16, fetch running"
    elif [ "$n_fetch" -ge "$MAX_FETCH" ]; then
      echo "$(date '+%F %T')  reads $n_ok/16 after $MAX_FETCH fetch submissions; giving up" >&2; exit 1
    else
      n_fetch=$((n_fetch+1))
      jid=$(sbatch --parsable scripts_oneoffs/malathion/fetch_sra.sh)
      echo "$(date '+%F %T')  reads $n_ok/16, submitted fetch $jid (#$n_fetch)"
    fi
    sleep 300; continue
  fi
  if [ ! -f "$CATALOG/catalog.tsv.gz" ]; then
    if squeue -u "$USER" -h -o %j | grep -q '^catalog_'; then
      echo "$(date '+%F %T')  reads 16/16, catalog building"
    else
      jid=$(bash pipeline/scripts/build_catalog.sh --founders pipeline/helpfiles/A_founders.bams.txt --out "$CATALOG" | tail -1)
      echo "$(date '+%F %T')  reads 16/16, submitted catalog build (last job $jid)"
    fi
    sleep 300; continue
  fi
  break
done

echo "$(date '+%F %T')  reads 16/16 and catalog present: running the standard pipeline"
bash scripts_oneoffs/malathion/run_2022.sh
