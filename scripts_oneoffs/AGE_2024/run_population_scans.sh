#!/bin/bash
# run_population_scans.sh — six female scans, split by source cage.
#
# Run from the repo root on HPC3, after `git pull`:
#   bash scripts_oneoffs/AGE_2024/run_population_scans.sh
#
# Both experiments drew from two cages built from the same eight founders, set up
# May 2023 and November 2023. The existing scans cannot separate them, and worse,
# AGE_SY replicates 1-6 -- the cut drawn in FigureX -- are entirely Nov while the
# pilot is 4/6 May. So food and population are confounded in that figure.
#
# Haplotypes already exist in both project folders, so this is scan-only: it
# submits smooth -> Wald+h2 -> concat for each design and nothing else.
#
#   lab_May    4 reps      lab_Nov     2 reps      -> process/AGE_2024
#   SY10_F_May 2           SY10_F_Nov 10           -> process/AGE_SY
#   SY20_F_May 2           SY20_F_Nov 10           -> process/AGE_SY
#
# Afterwards, re-run the gather (it picks these up) and fetch one file:
#   Rscript scripts_oneoffs/AGE_2024/gather_scans.R

set -euo pipefail

SMOOTH=100                # kb -- matches every other scan in this comparison
PILOT=process/AGE_2024
SY=process/AGE_SY

jobid_from() {
  local raw="$1" id
  id=$(printf '%s\n' "$raw" | tail -1 | tr -d '[:space:]')
  [[ "$id" =~ ^[0-9]+$ ]] || { echo "ERROR: no job id in:" >&2
                               printf '%s\n' "$raw" >&2; exit 1; }
  printf '%s' "$id"
}

[ -d "$PILOT/Haps" ] || { echo "ERROR: no haplotypes at $PILOT/Haps" >&2; exit 1; }
[ -d "$SY/Haps" ]    || { echo "ERROR: no haplotypes at $SY/Haps" >&2; exit 1; }

echo "building designs ..."
module load R/4.2.2 2>/dev/null || true
Rscript scripts_oneoffs/AGE_2024/make_population_designs.R

echo
echo "submitting ..."
submit() {   # name, dir, design
  local out; out=$(bash pipeline/scripts/run_scan.sh \
      --dir "$2" --scan "$1" --design "$3" --smooth "$SMOOTH")
  echo "   $(printf '%-18s' "$1") concat $(jobid_from "$out")"
}

submit AGE_2024_May "$PILOT" helpfiles/AGE_Aug13_24/Ageing_Aug13.May.txt
submit AGE_2024_Nov "$PILOT" helpfiles/AGE_Aug13_24/Ageing_Aug13.Nov.txt
for s in AGE_SY10_F AGE_SY20_F; do
  for p in May Nov; do
    submit "${s}_${p}" "$SY" "helpfiles/AGE_SY/${s}_${p}.test.txt"
  done
done

cat <<EOF

------------------------------------------------------------------
six scans submitted. when they finish:

  Rscript scripts_oneoffs/AGE_2024/gather_scans.R

then scp process/AGE_2024/AGE_2024_vs_AGE_SY.txt.gz as before.
------------------------------------------------------------------
EOF
