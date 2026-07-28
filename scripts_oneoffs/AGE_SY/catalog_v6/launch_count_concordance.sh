#!/bin/bash
# Gate step: runs AS A JOB after the catalog build finishes (call_samples self-submits,
# so it can't be a plain sbatch dependency). Counts the samples at the loosest catalog,
# then chains the concordance grid after the count merge.
#   args: <DIR> <CATDIR> <BAMLIST> <V3DIR>
set -euo pipefail
DIR=$1; CATDIR=$2; BAMLIST=$3; V3DIR=$4
JID=$(bash pipeline/scripts/call_samples.sh --catalog "$CATDIR" --bamlist "$BAMLIST" --dir "$DIR" | tail -1)
echo "count merge job: $JID"
[[ "$JID" =~ ^[0-9]+$ ]] || { echo "ERROR: no job id from call_samples (got '$JID')" >&2; exit 1; }
sbatch --dependency=afterok:"$JID" -o logs/AGE_SY/concordance_grid.out \
  scripts_oneoffs/AGE_SY/catalog_v6/concordance_grid.sh \
  "$DIR/Calls" "$V3DIR" "$CATDIR/snp_pass.tsv.gz"
echo "concordance chained after $JID"
