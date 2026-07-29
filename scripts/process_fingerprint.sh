#!/bin/bash
# process_fingerprint.sh — say what each process/<prefix>* folder actually is: stage
# (raw calls / counts / haplotypes / scans, QUAL layout OR catalog Catalog+Calls layout),
# how many samples and which, and newest mtime. Read-only. Writes logs/fingerprint.<prefix>.txt
# (syncs to Claude via cluster_sync).
#
#   bash scripts/process_fingerprint.sh AGE_SY
#   bash scripts/process_fingerprint.sh ZINC
set -uo pipefail
pat="${1:?usage: process_fingerprint.sh <prefix>   e.g. AGE_SY, ZINC}"
OUT="logs/fingerprint.${pat}.txt"; mkdir -p logs
{
shopt -s nullglob
for d in process/${pat}*/; do
  d=${d%/}
  echo "################ $d ################"
  echo "  newest mtime : $(find "$d" -maxdepth 3 -type f -printf '%TY-%Tm-%Td %TH:%TM\n' 2>/dev/null | sort | tail -1)"
  echo "  calls.bcf    : $(ls "$d"/calls.*.bcf 2>/dev/null | wc -l)"
  nl=$(find "$d" -maxdepth 1 -name 'RefAlt.chr*.txt' -type l 2>/dev/null | wc -l)
  nr=$(find "$d" -maxdepth 1 -name 'RefAlt.chr*.txt' -type f 2>/dev/null | wc -l)
  echo "  RefAlt (top) : $((nl+nr))  (real=$nr symlink=$nl)"
  echo "  RefAlt/founder (fromsites): $(find "$d" -maxdepth 1 -name 'RefAlt.*.chr*.txt' 2>/dev/null | grep -cvE '/RefAlt\.chr')"
  echo "  catalog layout: Catalog/=$([ -d "$d/Catalog" ] && echo yes || echo no)  Calls/RefAlt=$(find "$d/Calls" -maxdepth 1 -name 'RefAlt.chr*.txt' 2>/dev/null | wc -l)"
  echo "  R.haps       : $(ls "$d"/R.haps.*.rds 2>/dev/null | wc -l)"
  echo "  scan dirs    : $(find "$d" -maxdepth 1 -type d \( -name '*_F*' -o -name '*_M*' \) -printf '%f ' 2>/dev/null)"
  echo "  other        : $(ls -p "$d" 2>/dev/null | grep -viE '^(calls|RefAlt|R\.haps|Catalog/|Calls/)' | grep -vE '_[FM]/?$' | tr '\n' ' ')"
  # samples from a whole-sample RefAlt (top-level or Calls/), resolving symlink
  ra=$(ls "$d"/RefAlt.chrX.txt "$d"/Calls/RefAlt.chrX.txt 2>/dev/null | head -1)
  if [ -n "$ra" ]; then
    samps=$(head -1 "$ra" 2>/dev/null | tr ' \t' '\n' | sed -n 's/^REF_//p')
    n=$(printf '%s\n' "$samps" | grep -c .)
    echo "  samples ($n): $(printf '%s\n' "$samps" | head -8 | tr '\n' ' ')..."
  fi
  echo
done
} | tee "$OUT"
echo "wrote $OUT  ->  bash scripts/cluster_sync.sh"
