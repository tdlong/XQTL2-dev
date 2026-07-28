#!/bin/bash
# Gate step (runs as a job after catalog_filter): count the samples through the STANDARD
# pipeline (catalog is multiallelic-free now, so catalog_merge works), then chain the
# v3-vs-v6 comparison. No hack.
#   args: <DIR> <CATDIR> <BAMLIST> <V3DIR>
set -euo pipefail
DIR=$1; CATDIR=$2; BAMLIST=$3; V3DIR=$4
JID=$(bash pipeline/scripts/call_samples.sh --catalog "$CATDIR" --bamlist "$BAMLIST" --dir "$DIR" | tail -1)
echo "count merge job: $JID"
[[ "$JID" =~ ^[0-9]+$ ]] || { echo "ERROR: no job id from call_samples (got '$JID')" >&2; exit 1; }
sbatch --dependency=afterok:"$JID" -A tdlong_lab -p standard \
  --cpus-per-task=4 --mem-per-cpu=6G --time=00:40:00 -o logs/AGE_SY/compare_v3_v6.out \
  --wrap="module load R/4.2.2; Rscript pipeline/scripts/compare_refalt_calls.R $V3DIR $DIR/Calls"
echo "compare chained after $JID -> logs/AGE_SY/compare_v3_v6.out"
