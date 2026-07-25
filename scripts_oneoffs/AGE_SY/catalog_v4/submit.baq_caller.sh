#!/bin/bash
# One-command submit: array (4 settings, parallel) + dependent merge.
# Run once:  bash scripts_oneoffs/AGE_SY/catalog_v4/submit.baq_caller.sh
set -euo pipefail
S=scripts_oneoffs/AGE_SY/catalog_v4/baq_caller.chrX.sh
JID=$(sbatch --parsable --array=1-4 "$S")
MID=$(sbatch --parsable --dependency=afterok:"$JID" -o logs/AGE_SY/baq_caller_merge.out "$S" --merge)
echo "submitted: array $JID (tasks 1-4)  +  merge $MID (afterok:$JID)"
squeue -u "$USER"
