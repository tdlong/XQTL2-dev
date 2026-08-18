#!/bin/bash
# resume_haps_scan.sh — pick up AGE_2024 after the SNP calling, and finish it.
#
# The first run_all.sh submitted the 12 sample SNP calls and their merge fine,
# then failed to submit the haplotype job: call_samples.sh prints a banner line
# before the job ID, run_all.sh captured both lines, and the resulting
# --dependency was malformed, so sbatch rejected the haplotype job. Everything
# chained behind it sat as DependencyNeverSatisfied. run_all.sh is fixed; this
# script exists so the SNP calling does not have to be redone.
#
# The four AGE_SY reps 1-6 scans were unaffected and are running -- leave them.
#
# FIRST, cancel the three dead AGE_2024 jobs (they can never start):
#   scancel 55321193 55321194 55321195
#
# THEN, once the merge has finished (RefAlt.<chr>.txt exist), run:
#   bash scripts_oneoffs/AGE_2024/resume_haps_scan.sh
#
# It checks the merge output is really there rather than trusting the queue, so
# it is safe to run early -- it will just tell you to wait, and exit 75.
#
# To stop watching the queue altogether, leave this polling every 5 minutes and
# it will submit the moment the SNP calling lands:
#   until bash scripts_oneoffs/AGE_2024/resume_haps_scan.sh; do sleep 300; done
# (only "not ready yet" exits 75; a real error exits 1 and breaks the loop)

set -euo pipefail

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

# ── is the SNP calling actually finished? ────────────────────────────────────
missing=()
for c in "${CHRS[@]}"; do
  [ -s "$DIR/Calls/RefAlt.$c.txt" ] || missing+=("$c")
done
if [ ${#missing[@]} -gt 0 ]; then
  echo "The merge has not finished -- no RefAlt for: ${missing[*]}"
  echo "Counts present: $(ls "$DIR"/Calls/counts/*.tsv.gz 2>/dev/null | wc -l) of 20 (12 samples + 8 founders)"
  echo "Wait for the SNP calling (job 55321190) and its merge (55321191)."
  exit 75
fi

echo "RefAlt tables present for all 5 chromosomes:"
for c in "${CHRS[@]}"; do
  printf "   %-6s %s\n" "$c" "$(du -h "$DIR/Calls/RefAlt.$c.txt" | cut -f1)"
done

# sample columns should number 12 + 8 founders
ncol=$(head -1 "$DIR/Calls/RefAlt.chrX.txt" | tr '\t' '\n' | wc -l)
echo "   chrX header has $ncol fields"

# ── haplotypes, then the scan ────────────────────────────────────────────────
JID_HAPS=$(jobid_from "$(bash pipeline/scripts/run_haps.sh \
    --parfile "$PARFILE" --dir "$DIR")")
echo "haplotypes (75 kb / 5 kb) job: $JID_HAPS"

out=$(bash pipeline/scripts/run_scan.sh --after "$JID_HAPS" \
    --dir "$DIR" --scan AGE_2024 --design "$DESIGN" --smooth "$SMOOTH")
echo "AGE_2024 scan (smooth ${SMOOTH}kb):"
printf '%s\n' "$out" | sed 's/^/   /'

echo
echo "-> $DIR/Scans/AGE_2024/AGE_2024.scan.txt"
echo "then: Rscript scripts_oneoffs/AGE_2024/gather_scans.R"
