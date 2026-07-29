#!/bin/bash
# zinc_fingerprint.sh — say what each ZINC process folder actually is: dataset
# (A-lines/Hanson vs B-lines/ZINC2, from the RefAlt sample names), stage (calls / haps /
# scans / per-founder "fromsites" / downsample), sample count, newest mtime. Read-only.
# Writes logs/ZINC/fingerprint.txt (syncs to Claude via cluster_sync).
#   bash scripts_oneoffs/ZINC2/zinc_fingerprint.sh
set -uo pipefail
OUT=logs/ZINC/fingerprint.txt; mkdir -p logs/ZINC
DIRS="process/ZINC2 process/ZINC2_v2 process/ZINC2_fromsites process/ZINC2_fromsites_test process/ZINC_Hanson process/ZINC_Hanson_v2"
{
for d in $DIRS; do
  echo "################ $d ################"
  [ -d "$d" ] || { echo "  (absent)"; echo; continue; }
  echo "  newest mtime : $(find "$d" -maxdepth 2 -type f -printf '%TY-%Tm-%Td %TH:%TM\n' 2>/dev/null | sort | tail -1)"
  echo "  calls.bcf    : $(ls "$d"/calls.*.bcf 2>/dev/null | wc -l)"
  nl=$(find "$d" -maxdepth 1 -name 'RefAlt.chr*.txt' -type l 2>/dev/null | wc -l)
  nr=$(find "$d" -maxdepth 1 -name 'RefAlt.chr*.txt' -type f 2>/dev/null | wc -l)
  echo "  RefAlt whole : $((nl+nr))  (real=$nr symlink=$nl)"
  echo "  RefAlt /founder (fromsites): $(find "$d" -maxdepth 1 -name 'RefAlt.*.chr*.txt' 2>/dev/null | grep -cvE '/RefAlt\.chr')"
  echo "  R.haps       : $(ls "$d"/R.haps.*.rds 2>/dev/null | wc -l)"
  echo "  scan dirs    : $(find "$d" -maxdepth 1 -type d \( -name '*_F*' -o -name '*_M*' \) -printf '%f ' 2>/dev/null)"
  echo "  other files  : $(ls -p "$d" 2>/dev/null | grep -viE '^(calls|RefAlt|R\.haps)' | grep -vE '_[FM]/?$' | tr '\n' ' ')"
  ra="$d/RefAlt.chrX.txt"
  if [ -e "$ra" ]; then
    samps=$(head -1 "$ra" 2>/dev/null | tr ' \t' '\n' | sed -n 's/^REF_//p')
    n=$(printf '%s\n' "$samps" | grep -c .)
    echo "  samples ($n) : $(printf '%s\n' "$samps" | head -6 | tr '\n' ' ')..."
    case "$(printf '%s\n' "$samps" | head -1)" in
      C_*)  echo "  => A-lines / Hanson (samples C_x_y)";;
      Rep*) echo "  => B-lines / ZINC2 (samples RepNN_W/Z_M/F)";;
      *)    echo "  => (unrecognized sample naming)";;
    esac
  fi
  echo
done
} | tee "$OUT"
echo "wrote $OUT  ->  bash scripts/cluster_sync.sh"
