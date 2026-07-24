#!/usr/bin/env bash
# run_catalog_v4.sh — CALL AGE_SY samples against the founder catalog, then verify.
#
# New XQTL2 two-command design (main @ d40e279): building the catalog is a
# separate, deliberate one-time act; this script is the CALL step only.
#
# BUILD THE CATALOG FIRST (once per population; ~1.5 h; OVERWRITES --out):
#   bash pipeline/scripts/build_catalog.sh \
#       --founders pipeline/helpfiles/B_founders.bams.txt \
#       --out      process/AGE_SY_v4/Catalog
#
# THEN run this to count the samples (+ founders as columns) and auto-verify:
#   bash scripts_oneoffs/AGE_SY/catalog_v4/run_catalog_v4.sh
#
# To ADD samples later, point BAMLIST at a list of JUST the new BAMs and rerun —
# call_samples has no skip guard; new counts land next to the existing ones and
# RefAlt is re-merged. (Founders are already counted, so don't re-list them.)
#
# call_samples.sh submits its own SLURM array + merge and prints the merge JID; we
# chain the verification (compare_and_diagnose.sh) after it. Run from repo root.
set -uo pipefail

CATDIR=process/AGE_SY_v4/Catalog
DIR=process/AGE_SY_v4
BAMLIST=${BAMLIST:-helpfiles/AGE_SY/bam_list.v4.txt}   # 36 samples (R1-R6) + 8 founders

for f in "$BAMLIST" pipeline/scripts/call_samples.sh; do
  [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done
# The catalog must be built first — this step never builds it.
if [[ ! -s "$CATDIR/catalog.tsv.gz" ]]; then
  echo "ERROR: no catalog at $CATDIR/catalog.tsv.gz — build it first:" >&2
  echo "  bash pipeline/scripts/build_catalog.sh --founders pipeline/helpfiles/B_founders.bams.txt --out $CATDIR" >&2
  exit 1
fi

echo "calling $(grep -cve '^[[:space:]]*$' "$BAMLIST") BAMs from $BAMLIST against $CATDIR -> $DIR/Calls"
JID=$(bash pipeline/scripts/call_samples.sh \
        --catalog "$CATDIR" \
        --bamlist "$BAMLIST" \
        --dir     "$DIR")
echo "call_samples submitted; merge job id: $JID"

# chain the verification after the merge succeeds
SELFDIR=$(cd "$(dirname "$0")" && pwd)
VJID=$(sbatch --parsable --dependency=afterok:"$JID" "$SELFDIR/compare_and_diagnose.sh")
echo "verification chained after merge: job $VJID -> logs/AGE_SY/compare_v3_v4_${VJID}.out"
echo
echo "deliverable: $DIR/Calls/RefAlt.<chr>.txt"
echo "when everything finishes:  bash scripts/cluster_sync.sh"
