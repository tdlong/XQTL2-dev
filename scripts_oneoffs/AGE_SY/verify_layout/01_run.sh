#!/bin/bash
###############################################################################
# 01_run.sh — verify the XQTL2 #28 layout rewrite reproduces AGE_SY exactly.
#
# XQTL2 57042df moved pipeline output into Catalog/Calls/Haps/Scans and said:
#   "NEEDS VERIFICATION: I cannot run haps->smooth->scan here. XQTL2-dev to
#    confirm on AGE_SY that the scan still reproduces."
# This is that check. Nobody has run the new layout end to end.
#
# WHY THIS IS A CLEAN TEST
# The rewrite changed file paths and nothing else — no model, no parameter, no
# filter. So re-running the same inputs through the same code must reproduce the
# existing output BYTE FOR BYTE. Pass/fail, not judgement. 02_compare.sh checks it.
#
# WHAT IT DOES NOT TOUCH
# process/AGE_SY is the validated result and the baseline we compare against, so
# it is left completely alone. The test runs in its own folder that symlinks
# AGE_SY's Catalog/ and Calls/ read-only. Nothing is reorganized, moved, or
# overwritten. If this run is garbage, delete the folder and nothing is lost.
#
# Reproduces the run recorded in catalog_eval/NOTES.md:
#   REFALT2haps (size75k params) -> smooth 100 -> scan, four designs.
#
# Run from the repo root on the cluster:
#   bash scripts_oneoffs/AGE_SY/verify_layout/01_run.sh
# then when everything finishes:
#   bash scripts_oneoffs/AGE_SY/verify_layout/02_compare.sh
###############################################################################
set -euo pipefail

BASE=process/AGE_SY                       # validated result — READ ONLY
TEST=process/AGE_SY_layouttest            # this run's output
PARFILE=helpfiles/AGE_SY/AGE_SY_haplotype_parameters_size75k.R
SMOOTH=100
SCANS="AGE_SY10_F AGE_SY10_M AGE_SY20_F AGE_SY20_M"

# ── Preconditions ────────────────────────────────────────────────────────────
[ -d "$BASE/Calls" ]   || { echo "ERROR: $BASE/Calls not found" >&2; exit 1; }
[ -d "$BASE/Catalog" ] || { echo "ERROR: $BASE/Catalog not found" >&2; exit 1; }
[ -f "$PARFILE" ]      || { echo "ERROR: $PARFILE not found" >&2; exit 1; }
for c in chrX chr2L chr2R chr3L chr3R; do
  [ -f "$BASE/Calls/RefAlt.$c.txt" ]  || { echo "ERROR: missing $BASE/Calls/RefAlt.$c.txt" >&2; exit 1; }
  [ -f "$BASE/R.haps.$c.out.rds" ]    || { echo "ERROR: missing baseline haps $BASE/R.haps.$c.out.rds" >&2; exit 1; }
done
for s in $SCANS; do
  [ -f "helpfiles/AGE_SY/$s.test.txt" ]  || { echo "ERROR: design missing: helpfiles/AGE_SY/$s.test.txt" >&2; exit 1; }
  [ -f "$BASE/$s/$s.scan.txt" ]          || { echo "ERROR: baseline scan missing: $BASE/$s/$s.scan.txt" >&2; exit 1; }
done
[ -e "$TEST" ] && { echo "ERROR: $TEST already exists — remove it first (rm -rf $TEST)" >&2; exit 1; }

# ── Test folder: Catalog/ and Calls/ symlinked from the real project ─────────
# Absolute targets so the links survive regardless of where they are read from.
mkdir -p "$TEST"
ln -s "$(cd "$BASE/Catalog" && pwd)" "$TEST/Catalog"
ln -s "$(cd "$BASE/Calls"   && pwd)" "$TEST/Calls"
echo "test folder: $TEST  (Catalog/ and Calls/ -> $BASE, read only)"

# ── Step 4: haplotypes -> $TEST/Haps/ ────────────────────────────────────────
# run_haps.sh already asks for highmem 10G / 1 day / array 1-5 (Slurm.md: highmem
# is 10 GB per core, so --mem-per-cpu=10G is the max without pulling extra cores).
JID_HAPS=$(bash pipeline/scripts/run_haps.sh --parfile "$PARFILE" --dir "$TEST")
echo "haplotypes submitted: $JID_HAPS  -> $TEST/Haps/"

# ── Step 5: the four scans -> $TEST/Scans/<scan>/ ────────────────────────────
# Submitted in parallel, all chained behind the haplotype array, so the extra
# three scans cost disk but no wall time.
for s in $SCANS; do
  out=$(bash pipeline/scripts/run_scan.sh \
        --after   "$JID_HAPS" \
        --design  "helpfiles/AGE_SY/$s.test.txt" \
        --dir     "$TEST" \
        --scan    "$s" \
        --smooth  "$SMOOTH")
  echo "$out" | sed "s/^/  [$s] /"
done

cat <<EOF

================================================================
Submitted as one chain. Nothing else to run.
  haps ($JID_HAPS) -> smooth ${SMOOTH}kb -> r2 -> scan -> concat, x4 designs

$BASE is untouched — it is the baseline.
Output lands in $TEST/{Haps,Scans}/

Monitor:  squeue -u \$USER
When everything has finished:
  bash scripts_oneoffs/AGE_SY/verify_layout/02_compare.sh
================================================================
EOF
