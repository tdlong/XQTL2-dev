#!/bin/bash
###############################################################################
# 01_rebuild_AGE_SY.sh — rebuild AGE_SY under the #28 layout, keeping the old
# one alongside so we can prove the rewrite reproduces it.
#
# WHY
# XQTL2 57042df moved pipeline output into Catalog/Calls/Haps/Scans and left a
# note: "NEEDS VERIFICATION: I cannot run haps->smooth->scan here. XQTL2-dev to
# confirm on AGE_SY that the scan still reproduces." Nobody has run it. AGE_SY is
# the only project with finished results to check against.
#
# WHAT THIS DOES (all on the cluster)
#   1. process/AGE_SY -> process/AGE_SY_flatorg          the old result, parked
#   2. a fresh, empty process/AGE_SY
#   3. copy Catalog/ and Calls/ (RefAlt + counts) into it
#   4. run step 4 (haplotypes) and step 5a (smooth -> scan) on the new AGE_SY
#
# Steps 1-3 of the pipeline are NOT repeated: the catalog and the per-sample
# counts are copied, not rebuilt. Only the two steps whose code #28 changed get
# re-run.
#
# The copy deliberately skips Calls/compare.* — 2.07 GB of leftover debug output
# from the v3-vs-v6 comparison. Nothing reads it and it is not pipeline output.
# It stays in AGE_SY_flatorg and goes away when that folder does.
#
# Parameters reproduce the run recorded in catalog_eval/NOTES.md:
# size75k haplotype params, smooth 100, the four designs.
#
# Run from the repo root on the cluster:
#   bash scripts_oneoffs/AGE_SY/verify_layout/01_rebuild_AGE_SY.sh
# then when the jobs finish:
#   bash scripts_oneoffs/AGE_SY/verify_layout/02_compare.sh
#
# Nothing is deleted. If this goes wrong, undo is:
#   rm -rf process/AGE_SY && mv process/AGE_SY_flatorg process/AGE_SY
###############################################################################
set -euo pipefail

NEW=process/AGE_SY
OLD=process/AGE_SY_flatorg
PARFILE=helpfiles/AGE_SY/AGE_SY_haplotype_parameters_size75k.R
SMOOTH=100
CHRS="chrX chr2L chr2R chr3L chr3R"
SCANS="AGE_SY10_F AGE_SY10_M AGE_SY20_F AGE_SY20_M"

# ── Preconditions — check everything BEFORE renaming anything ────────────────
[ -e "$OLD" ] && { echo "ERROR: $OLD already exists. Previous attempt? Resolve it first." >&2; exit 1; }
[ -d "$NEW/Calls" ]   || { echo "ERROR: $NEW/Calls not found — is this the right repo root?" >&2; exit 1; }
[ -d "$NEW/Catalog" ] || { echo "ERROR: $NEW/Catalog not found" >&2; exit 1; }
[ -f "$PARFILE" ]     || { echo "ERROR: $PARFILE not found" >&2; exit 1; }
for c in $CHRS; do
  [ -f "$NEW/Calls/RefAlt.$c.txt" ] || { echo "ERROR: missing $NEW/Calls/RefAlt.$c.txt" >&2; exit 1; }
  [ -f "$NEW/R.haps.$c.out.rds" ]   || { echo "ERROR: missing baseline haps $NEW/R.haps.$c.out.rds" >&2; exit 1; }
done
for s in $SCANS; do
  [ -f "helpfiles/AGE_SY/$s.test.txt" ] || { echo "ERROR: design missing: helpfiles/AGE_SY/$s.test.txt" >&2; exit 1; }
  [ -f "$NEW/$s/$s.scan.txt" ]          || { echo "ERROR: baseline scan missing: $NEW/$s/$s.scan.txt" >&2; exit 1; }
done
[ -d "$NEW/Calls/counts" ] || { echo "ERROR: $NEW/Calls/counts not found" >&2; exit 1; }
echo "preconditions OK — baseline haps and all four baseline scans are present."

# ── 1. park the old result ───────────────────────────────────────────────────
mv "$NEW" "$OLD"
echo "parked: $NEW -> $OLD"

# ── 2-3. fresh AGE_SY with the calling stages copied in ──────────────────────
# cp -a preserves timestamps, so the copied catalog/counts still show when they
# were actually made rather than when they were copied.
mkdir -p "$NEW/Calls"
cp -a "$OLD/Catalog"              "$NEW/Catalog"
cp -a "$OLD/Calls/counts"         "$NEW/Calls/counts"
cp -a "$OLD"/Calls/RefAlt.*.txt   "$NEW/Calls/"
echo "copied Catalog/ and Calls/ (RefAlt + counts; compare.* debug files skipped)"

# verify the copy before spending cluster time on it
for c in $CHRS; do
  a=$(stat -c%s "$OLD/Calls/RefAlt.$c.txt"); b=$(stat -c%s "$NEW/Calls/RefAlt.$c.txt")
  [ "$a" = "$b" ] || { echo "ERROR: RefAlt.$c.txt copied short ($a vs $b)" >&2; exit 1; }
done
na=$(ls "$OLD/Calls/counts" | wc -l); nb=$(ls "$NEW/Calls/counts" | wc -l)
[ "$na" = "$nb" ] || { echo "ERROR: counts/ copied $nb of $na files" >&2; exit 1; }
echo "copy verified: 5 RefAlt tables, $nb count files, catalog."

# ── 4. step 4 — haplotypes -> $NEW/Haps/ ─────────────────────────────────────
# run_haps.sh asks for highmem / 10G per core / array 1-5 on its own (Slurm.md:
# highmem is 10 GB per core, the max without implicitly pulling extra cores).
JID_HAPS=$(bash pipeline/scripts/run_haps.sh --parfile "$PARFILE" --dir "$NEW")
echo "haplotypes submitted: $JID_HAPS  -> $NEW/Haps/"

# ── 5a. the four scans -> $NEW/Scans/<scan>/ ─────────────────────────────────
# All four chained behind the haplotype array and run in parallel with each other.
for s in $SCANS; do
  bash pipeline/scripts/run_scan.sh \
      --after   "$JID_HAPS" \
      --design  "helpfiles/AGE_SY/$s.test.txt" \
      --dir     "$NEW" \
      --scan    "$s" \
      --smooth  "$SMOOTH" | sed "s/^/  [$s] /"
done

cat <<EOF

================================================================
Submitted. Nothing else to run.
  haps ($JID_HAPS) -> smooth ${SMOOTH}kb -> r2 -> scan -> concat, x4 designs

  $NEW   being rebuilt under the #28 layout
  $OLD   the old result, untouched — this is the answer sheet

Monitor:  squeue -u \$USER
When everything has finished:
  bash scripts_oneoffs/AGE_SY/verify_layout/02_compare.sh
  bash scripts/cluster_sync.sh

Delete $OLD only after 02_compare.sh says PASS.
================================================================
EOF
