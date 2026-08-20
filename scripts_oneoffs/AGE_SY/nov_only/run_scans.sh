#!/bin/bash
# run_scans.sh — the twelve AGE_SY scans without replicates 8 and 9.
#
# Run from the repo root on HPC3, after `git pull`:
#   bash scripts_oneoffs/AGE_SY/nov_only/run_scans.sh
#
# Replicates 8 and 9 are the May 2023 cage; dropping them leaves a single-cage
# experiment on 10 replicates, splitting evenly 5/5 into odd and even.
#
#   4 x <scan>_no89              -> process/AGE_SY/Scans            (10 reps)
#   8 x <scan>_no89_{odd,even}   -> process/AGE_SY_splithalf/Scans  (5 reps each)
#
# Haplotypes already exist in both folders, so this is scan-only. The split-half
# folder reads its Haps through a symlink to process/AGE_SY, set up by
# scripts_oneoffs/AGE_SY/splithalf/run_splithalf_scans.sh -- it must already be
# there, and this script checks.
#
# Afterwards, on the cluster:
#   Rscript scripts_oneoffs/AGE_SY/nov_only/gather.R
# then fetch the two files it names.

set -euo pipefail

SMOOTH=100                       # kb -- matches every other AGE_SY scan
FULL_DIR=process/AGE_SY
HALF_DIR=process/AGE_SY_splithalf
DESIGNS=helpfiles/AGE_SY/nov_only
SCANS=(AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M)

jobid_from() {
  local id
  id=$(printf '%s\n' "$1" | grep -oE '[0-9]+' | tail -1)
  [[ -n "$id" ]] || { echo "ERROR: no job id in:" >&2
                      printf '%s\n' "$1" >&2; return 1; }
  printf '%s' "$id"
}

[ -d "$FULL_DIR/Haps" ] || { echo "ERROR: no haplotypes at $FULL_DIR/Haps" >&2; exit 1; }
[ -e "$HALF_DIR/Haps" ] || { echo "ERROR: no Haps symlink at $HALF_DIR/Haps" >&2; exit 1; }

echo "building designs ..."
module load R/4.2.2 2>/dev/null || true
Rscript scripts_oneoffs/AGE_SY/nov_only/make_designs.R

echo
echo "submitting ..."
submit() {   # name, dir, design
  [ -f "$3" ] || { echo "ERROR: missing design $3" >&2; exit 1; }
  local out; out=$(bash pipeline/scripts/run_scan.sh \
      --dir "$2" --scan "$1" --design "$3" --smooth "$SMOOTH")
  echo "   $(printf '%-26s' "$1") concat $(jobid_from "$out")"
}

for s in "${SCANS[@]}"; do
  submit "${s}_no89" "$FULL_DIR" "$DESIGNS/${s}.no89.txt"
done
for s in "${SCANS[@]}"; do
  for h in odd even; do
    submit "${s}_no89_${h}" "$HALF_DIR" "$DESIGNS/${s}.no89.${h}.txt"
  done
done

cat <<EOF

------------------------------------------------------------------
12 scans submitted (4 x 10 reps, 8 x 5 reps). when they finish:

  Rscript scripts_oneoffs/AGE_SY/nov_only/gather.R

which writes two files to fetch:
  $FULL_DIR/AGE_SY_4scan_no89.txt.gz
  $HALF_DIR/AGE_SY_splithalf_H2_no89.txt.gz
------------------------------------------------------------------
EOF
