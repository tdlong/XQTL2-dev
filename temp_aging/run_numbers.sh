#!/bin/bash
# run_numbers.sh — run every analysis script and record its output WITH the
# provenance of the data that produced it.
#
#   bash temp_aging/run_numbers.sh
#
# Writes temp_aging/numbers/<script>.txt. Each file starts with a header naming
# every input the script read, that input's size, modification time and MD5, and
# the git commit of the script itself. So any number quoted in the prose can be
# traced to the exact bytes it came from.
#
# The reason this exists: the scripts printed to stdout and nothing was saved, so
# a number in METHODS or a legend had no link back to a dataset. There was no way
# to tell a 12-replicate number from a 10-replicate one by looking at it.

set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

OUT=temp_aging/numbers
mkdir -p "$OUT"

SCRIPTS=(h2_threshold scan_resolution significant_regions peak_table chr3L_peak)

stamp_inputs() {   # script path -> lines describing each process/ file it reads
  grep -ohE '"process/[^"]+\.txt\.gz"' "$1" | tr -d '"' | sort -u | while read -r f; do
    if [ -f "$f" ]; then
      printf '#   %-52s %8s  %s  md5:%s\n' "$f" \
        "$(du -h "$f" | cut -f1)" \
        "$(date -r "$f" '+%Y-%m-%d %H:%M')" \
        "$(md5 -q "$f" 2>/dev/null | cut -c1-8 || md5sum "$f" | cut -c1-8)"
    else
      printf '#   %-52s MISSING\n' "$f"
    fi
  done
}

for s in "${SCRIPTS[@]}"; do
  src="temp_aging/$s.R"
  dst="$OUT/$s.txt"
  [ -f "$src" ] || { echo "missing $src" >&2; continue; }
  {
    echo "# $s.R"
    echo "# run:    $(date '+%Y-%m-%d %H:%M')"
    echo "# script: $(git log -1 --format=%h -- "$src" 2>/dev/null || echo uncommitted)"
    echo "# reads:"
    stamp_inputs "$src"
    echo "#"
  } > "$dst"
  Rscript "$src" 2>&1 | grep -vE '^(Warning message|In addition|package .* was built|  [0-9]+:)' >> "$dst"
  printf "  %-26s -> %s\n" "$s.R" "$dst"
done

echo
echo "every number in the prose should be findable in one of these files."
