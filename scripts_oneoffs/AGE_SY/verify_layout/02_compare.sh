#!/bin/bash
###############################################################################
# 02_compare.sh — did the new layout reproduce AGE_SY exactly?
#
# Run after 01_run.sh finishes. Compares the test run against the validated
# result. The rewrite changed only file paths, so everything must match exactly;
# any difference is a bug in the layout rewrite, not a scientific finding.
#
#   bash scripts_oneoffs/AGE_SY/verify_layout/02_compare.sh
#
# Writes logs/verify_layout.txt (syncs to Claude via cluster_sync).
###############################################################################
set -uo pipefail

BASE=process/AGE_SY
TEST=process/AGE_SY_layouttest
CHRS="chrX chr2L chr2R chr3L chr3R"
SCANS="AGE_SY10_F AGE_SY10_M AGE_SY20_F AGE_SY20_M"
OUT=logs/verify_layout.txt; mkdir -p logs

fail=0
{
echo "AGE_SY layout verification — $(date)"
echo "baseline: $BASE     test: $TEST"
echo

# ── The scans are the real test: plain text, fully deterministic ─────────────
echo "=== scan output (the authoritative check) ==="
for s in $SCANS; do
  a="$BASE/$s/$s.scan.txt"
  b="$TEST/Scans/$s/$s.scan.txt"
  if [ ! -f "$b" ]; then
    echo "  MISSING  $b   (scan did not finish — check squeue / logs)"; fail=1; continue
  fi
  if cmp -s "$a" "$b"; then
    echo "  IDENTICAL  $s.scan.txt"
  else
    echo "  DIFFERS    $s.scan.txt"; fail=1
    echo "             baseline $(wc -l < "$a") lines, test $(wc -l < "$b") lines"
    echo "             first differing line:"
    diff "$a" "$b" 2>/dev/null | head -4 | sed 's/^/               /'
  fi
done
echo

# ── Haplotypes: .rds is gzip-compressed, so a byte diff is not proof of a real
# ── difference. Report it, but let the scan comparison above be the verdict.
echo "=== haplotypes (informational — .rds compression need not be byte-stable) ==="
for c in $CHRS; do
  a="$BASE/R.haps.$c.out.rds"
  b="$TEST/Haps/R.haps.$c.out.rds"
  if [ ! -f "$b" ]; then
    echo "  MISSING  $b"; fail=1; continue
  fi
  sa=$(stat -c%s "$a" 2>/dev/null || stat -f%z "$a")
  sb=$(stat -c%s "$b" 2>/dev/null || stat -f%z "$b")
  if cmp -s "$a" "$b"; then echo "  IDENTICAL  R.haps.$c.out.rds"
  else                      echo "  byte-differs  R.haps.$c.out.rds  (baseline ${sa} B, test ${sb} B)"
  fi
done
echo

# ── Confirm the files actually landed where #28 says they should ─────────────
echo "=== layout check: did the pipeline write to the new locations? ==="
for d in Haps Scans; do
  [ -d "$TEST/$d" ] && echo "  OK       $TEST/$d/ exists" || { echo "  MISSING  $TEST/$d/"; fail=1; }
done
for s in $SCANS; do
  [ -d "$TEST/Scans/$s" ] && echo "  OK       $TEST/Scans/$s/" || { echo "  MISSING  $TEST/Scans/$s/"; fail=1; }
done
stray=$(find "$TEST" -maxdepth 1 -name 'R.haps.*' -o -maxdepth 1 -name 'RefAlt.*' 2>/dev/null)
[ -n "$stray" ] && { echo "  UNEXPECTED files at top level (should be under Haps/ or Calls/):"; \
                     printf '%s\n' "$stray" | sed 's/^/    /'; fail=1; }
echo

if [ "$fail" -eq 0 ]; then
  echo "RESULT: PASS — the #28 layout reproduces AGE_SY exactly."
  echo "        Safe to report back on XQTL2 57042df / #28."
else
  echo "RESULT: FAIL — see above. Do not migrate other projects until this is understood."
fi
} | tee "$OUT"
echo
echo "wrote $OUT  ->  bash scripts/cluster_sync.sh"
