#!/bin/bash
# submit.call_v5.sh — count the AGE_SY samples against the v5 catalog. One command:
#   bash scripts_oneoffs/AGE_SY/catalog_v5/submit.call_v5.sh
#
# Same 36-sample BAM list as v4 (helpfiles/AGE_SY/bam_list.v4.txt, R1-R6 + 8 founders)
# so v4-vs-v5 differ ONLY in the caller/catalog version. call_samples.sh submits its
# own SLURM array + merge and returns the merge job id. Deliverable:
#   process/AGE_SY_v5/Calls/RefAlt.<chr>.txt   (post-#14-fix -> clean v3-vs-v5 compare)
set -euo pipefail
CATDIR=process/AGE_SY_v5/Catalog
DIR=process/AGE_SY_v5
BAMLIST=${BAMLIST:-helpfiles/AGE_SY/bam_list.v4.txt}

[[ -s "$CATDIR/catalog.tsv.gz" ]] || { echo "ERROR: no v5 catalog at $CATDIR/catalog.tsv.gz — build it first." >&2; exit 1; }
[[ -s "$BAMLIST" ]] || { echo "ERROR: missing BAM list $BAMLIST" >&2; exit 1; }

echo "calling $(grep -cve '^[[:space:]]*$' "$BAMLIST") BAMs from $BAMLIST against $CATDIR -> $DIR/Calls"
JID=$(bash pipeline/scripts/call_samples.sh --catalog "$CATDIR" --bamlist "$BAMLIST" --dir "$DIR")
echo "v5 call submitted; merge job id: $JID  -> $DIR/Calls/RefAlt.<chr>.txt"
squeue -u "$USER"
