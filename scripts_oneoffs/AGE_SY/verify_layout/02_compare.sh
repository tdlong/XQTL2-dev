#!/bin/bash
###############################################################################
# 02_compare.sh — did the rebuilt AGE_SY reproduce the old one?
#
# Run after 01_rebuild_AGE_SY.sh finishes.
#
#   process/AGE_SY           rebuilt under the #28 layout   (Haps/, Scans/<scan>/)
#   process/AGE_SY_flatorg   the old result                 (R.haps.*, <scan>/)
#
# #28 changed file paths and nothing else — no model, no parameter, no filter.
# Same inputs, same parameters, so the output must match exactly. A difference is
# a bug in the rewrite, not a result.
#
#   bash scripts_oneoffs/AGE_SY/verify_layout/02_compare.sh
#
# Writes logs/verify_layout.txt, which cluster_sync.sh pushes back.
###############################################################################
set -uo pipefail

NEW=process/AGE_SY
OLD=process/AGE_SY_flatorg
CHRS="chrX chr2L chr2R chr3L chr3R"
SCANS="AGE_SY10_F AGE_SY10_M AGE_SY20_F AGE_SY20_M"
OUT=logs/verify_layout.txt; mkdir -p logs

fail=0
{
echo "AGE_SY #28 layout verification — $(date)"
echo "  new (rebuilt): $NEW"
echo "  old (baseline): $OLD"
echo

# ── The scans decide it: plain text, fully deterministic ─────────────────────
echo "=== scan output — this is the verdict ==="
for s in $SCANS; do
  a="$OLD/$s/$s.scan.txt"
  b="$NEW/Scans/$s/$s.scan.txt"
  if [ ! -f "$b" ]; then
    echo "  MISSING    $b   (scan did not finish — check squeue and logs/)"; fail=1; continue
  fi
  if cmp -s "$a" "$b"; then
    echo "  IDENTICAL  $s.scan.txt"
  else
    echo "  DIFFERS    $s.scan.txt"; fail=1
    echo "             old $(wc -l < "$a") lines, new $(wc -l < "$b") lines"
    diff "$a" "$b" 2>/dev/null | head -4 | sed 's/^/               /'
  fi
done
echo

# ── Haplotypes: .rds is gzip-compressed, so a byte difference is not by itself
# ── proof of a real difference. Reported, but the scans above are the verdict.
echo "=== haplotypes (informational — .rds need not be byte-stable) ==="
for c in $CHRS; do
  a="$OLD/R.haps.$c.out.rds"
  b="$NEW/Haps/R.haps.$c.out.rds"
  if [ ! -f "$b" ]; then echo "  MISSING    $b"; fail=1; continue; fi
  if cmp -s "$a" "$b"; then
    echo "  IDENTICAL  R.haps.$c.out.rds"
  else
    echo "  byte-differs  R.haps.$c.out.rds  (old $(stat -c%s "$a") B, new $(stat -c%s "$b") B)"
  fi
done
echo

# ── Did the pipeline actually write to the new locations? ────────────────────
echo "=== layout: are the files where #28 says they should be? ==="
for d in Catalog Calls Haps Scans; do
  [ -d "$NEW/$d" ] && echo "  OK         $NEW/$d/" || { echo "  MISSING    $NEW/$d/"; fail=1; }
done
for s in $SCANS; do
  [ -d "$NEW/Scans/$s" ] && echo "  OK         $NEW/Scans/$s/" || { echo "  MISSING    $NEW/Scans/$s/"; fail=1; }
done
stray=$(find "$NEW" -maxdepth 1 \( -name 'R.haps.*' -o -name 'RefAlt.*' \) 2>/dev/null)
if [ -n "$stray" ]; then
  echo "  UNEXPECTED at top level (should be under Haps/ or Calls/):"
  printf '%s\n' "$stray" | sed 's/^/    /'; fail=1
fi
echo

if [ "$fail" -eq 0 ]; then
  echo "RESULT: PASS — the #28 layout reproduces AGE_SY exactly."
  echo "        Safe to delete the baseline:  rm -rf $OLD"
  echo "        Report back on XQTL2 57042df / #28."
else
  echo "RESULT: FAIL — see above."
  echo "        KEEP $OLD. Do not migrate any other project until this is understood."
  echo "        Full undo:  rm -rf $NEW && mv $OLD $NEW"
fi
} | tee "$OUT"
echo
echo "wrote $OUT  ->  bash scripts/cluster_sync.sh"
