#!/usr/bin/env bash
# status.sh — where is this analysis, right now.
#
#   bash temp_aging/status.sh
#
# Run from the XQTL2-dev root, on HPC3 or the laptop. Read-only, seconds.
# Answers: which commit am I on, which products exist, how old are they
# relative to the scans they derive from, and is anything still queued.

set -uo pipefail
[ -d scripts_oneoffs/AGE_SY ] || { echo "run from the XQTL2-dev root" >&2; exit 1; }

age () { [ -f "$1" ] && date -r "$1" '+%Y-%m-%d %H:%M' 2>/dev/null || echo "MISSING"; }
row () { printf '  %-56s %s\n' "$1" "$(age "$2")"; }

echo "=== repo ==="
printf '  XQTL2-dev  %s  %s\n' "$(git log -1 --format=%h)" "$(git log -1 --format=%s | cut -c1-58)"
[ -d pipeline ] && printf '  pipeline   %s  %s\n' \
  "$(git -C pipeline log -1 --format=%h)" "$(git -C pipeline log -1 --format=%s | cut -c1-58)"
git fetch -q origin 2>/dev/null && \
  printf '  behind origin by %s commit(s)\n' "$(git rev-list --count HEAD..origin/main 2>/dev/null || echo '?')"

echo
echo "=== the scans (primary; everything else derives from these) ==="
for sc in AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M; do
  d=process/AGE_SY/Scans/${sc}_no89
  row "${sc}_no89 scan.txt"    "$d/${sc}_no89.scan.txt"
  row "${sc}_no89 smooth_r2"   "$d/${sc}_no89.smooth_r2.txt"
done

echo
echo "=== derived ==="
row "4scan"        process/AGE_SY/AGE_SY_4scan_no89.txt.gz
row "splithalf H2" process/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz
row "varcomp (Fig 1c)" process/AGE_SY_splithalf/H2_varcomp_by_window_no89.txt.gz
row "zoom means (Fig 2)" process/AGE_SY/AGE_SY_zoom_means.txt.gz
row "SNP scan (Fig 3)" process/AGE_SY/AGE_SY_4snpscan_no89.txt.gz
row "numbers/coverage" temp_aging/numbers/coverage.txt
row "peak_table"      temp_aging/peak_table.txt

echo
echo "=== stale? (derived older than the newest scan is a problem) ==="
newest=$(ls -t process/AGE_SY/Scans/*/*.scan.txt 2>/dev/null | head -1)
if [ -n "$newest" ]; then
  echo "  newest scan.txt: $(age "$newest")"
  for f in process/AGE_SY/AGE_SY_4scan_no89.txt.gz \
           process/AGE_SY_splithalf/H2_varcomp_by_window_no89.txt.gz \
           process/AGE_SY/AGE_SY_zoom_means.txt.gz \
           process/AGE_SY/AGE_SY_4snpscan_no89.txt.gz; do
    if [ ! -f "$f" ]; then printf '  MISSING  %s\n' "$f"
    elif [ "$f" -ot "$newest" ]; then printf '  STALE    %s\n' "$f"
    fi
  done
else
  echo "  no scan.txt found"
fi

echo
if command -v squeue >/dev/null 2>&1; then
  echo "=== queue ==="
  squeue -u "$USER" -o '%.10i %.20j %.2t %.10M %R' 2>/dev/null | head -20
fi
