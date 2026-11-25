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
# there runs the standard pipeline (run_full_pipeline.sh). Waits for that to
# finish, then runs the per-pair scans and the paper's ANOVA on the new
# haplotypes (rerun_2022.sh -> finish_2022.sh, which pushes the results to
# logs/malathion/), and finally syncs its own log. Nothing to do by hand.
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
out=$(bash scripts_oneoffs/malathion/run_2022.sh) || { echo "$out"; echo "run_2022.sh failed" >&2; bash scripts/cluster_sync.sh; exit 1; }
echo "$out"
final=$(printf '%s\n' "$out" | sed -n 's/.*Final job: \([0-9]*\).*/\1/p' | tail -1)
[ -n "$final" ] || { echo "could not find the final job id in run_full_pipeline output" >&2; bash scripts/cluster_sync.sh; exit 1; }

# wait for the standard pipeline (figures job is last)
while squeue -h -j "$final" 2>/dev/null | grep -q .; do
  echo "$(date '+%F %T')  pipeline running (final job $final)"; sleep 600
done
SCAN=process/malathion_2022/Scans/malathion_2022_all/malathion_2022_all.scan.txt
if [ ! -f "$SCAN" ] || [ ! -f process/malathion_2022/Haps/R.haps.chr3R.out.rds ]; then
  echo "$(date '+%F %T')  pipeline ended without $SCAN; see the slurm-*.out logs" >&2
  bash scripts/cluster_sync.sh; exit 1
fi
echo "$(date '+%F %T')  standard pipeline done: $SCAN"

# per-pair scans + the paper's ANOVA on these haplotypes; finish_2022.sh pushes
mkdir -p logs/malathion
cp process/malathion_2022/Scans/malathion_2022_all/*.png logs/malathion/ 2>/dev/null
bash scripts_oneoffs/malathion/rerun_2022.sh \
    --refalt process/malathion_2022/Calls \
    --reuse-haps process/malathion_2022/Haps \
    --design helpfiles/malathion/design_2022.txt \
    --parfile helpfiles/malathion/hap_params_2022.R \
    --dir process/malathion_2022_pairs \
    --prefix malathion_2022
echo "$(date '+%F %T')  diagnostics submitted; finish job pushes results. syncing this log."
bash scripts/cluster_sync.sh
