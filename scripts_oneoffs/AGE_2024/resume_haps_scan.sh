#!/bin/bash
# resume_haps_scan.sh — submit the two AGE_2024 steps that never got queued.
#
# The first run_all.sh submitted the 12 sample SNP calls and their merge fine,
# then failed to submit the haplotype job: call_samples.sh prints a banner line
# before the job ID, run_all.sh captured both lines, and the resulting
# --dependency was malformed so sbatch rejected it. Everything behind it sat as
# DependencyNeverSatisfied. run_all.sh is fixed; this finishes the run without
# redoing the SNP calling.
#
# Cancel the three dead jobs first -- they can never start:
#   scancel 55321193 55321194 55321195
#
# Then submit, chained on the merge job so there is nothing to wait for:
#   bash scripts_oneoffs/AGE_2024/resume_haps_scan.sh 55321191
#
# With no argument it requires the RefAlt tables to already be on disk.
#
# The four AGE_SY reps 1-6 scans were unaffected -- leave them alone.

set -euo pipefail

AFTER="${1:-}"
DIR=process/AGE_2024
PARFILE=helpfiles/AGE_Aug13_24/AGE_2024_haplotype_parameters_size75k.R
DESIGN=helpfiles/AGE_Aug13_24/Ageing_Aug13.txt
SMOOTH=100
CHRS=(chrX chr2L chr2R chr3L chr3R)

jobid_from() {
  local raw="$1" id
  id=$(printf '%s\n' "$raw" | tail -1 | tr -d '[:space:]')
  [[ "$id" =~ ^[0-9]+$ ]] || { echo "ERROR: no job id in:" >&2
                               printf '%s\n' "$raw" >&2; exit 1; }
  printf '%s' "$id"
}

for f in "$PARFILE" "$DESIGN"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

if [ -n "$AFTER" ]; then
  [[ "$AFTER" =~ ^[0-9]+$ ]] || { echo "ERROR: '$AFTER' is not a job id" >&2; exit 1; }
  AFTER_FLAG=(--after "$AFTER")
  echo "chaining on job $AFTER (the SNP-call merge)"
else
  missing=()
  for c in "${CHRS[@]}"; do
    [ -s "$DIR/Calls/RefAlt.$c.txt" ] || missing+=("$c")
  done
  if [ ${#missing[@]} -gt 0 ]; then
    echo "No RefAlt table for: ${missing[*]}" >&2
    echo "Either wait for the merge, or pass its job id:" >&2
    echo "  bash $0 <merge_jobid>" >&2
    exit 1
  fi
  AFTER_FLAG=()
  echo "RefAlt tables already present -- submitting immediately"
fi

JID_HAPS=$(jobid_from "$(bash pipeline/scripts/run_haps.sh "${AFTER_FLAG[@]}" \
    --parfile "$PARFILE" --dir "$DIR")")
echo "haplotypes (75 kb / 5 kb): $JID_HAPS"

out=$(bash pipeline/scripts/run_scan.sh --after "$JID_HAPS" \
    --dir "$DIR" --scan AGE_2024 --design "$DESIGN" --smooth "$SMOOTH")
printf '%s\n' "$out" | sed 's/^/   /'

echo
echo "-> $DIR/Scans/AGE_2024/AGE_2024.scan.txt"
echo "then: Rscript scripts_oneoffs/AGE_2024/gather_scans.R"
