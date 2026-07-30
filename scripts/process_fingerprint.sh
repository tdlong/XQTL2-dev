#!/bin/bash
# process_fingerprint.sh — say what each process/<prefix>* folder actually is: which
# layout it uses, what stages exist (calls / counts / haplotypes / scans), how many
# samples and which, and newest mtime. Read-only. Writes logs/fingerprint.<prefix>.txt
# (syncs to Claude via cluster_sync).
#
#   bash scripts/process_fingerprint.sh AGE_SY
#   bash scripts/process_fingerprint.sh ZINC
#
# Two layouts exist. XQTL2 #28 moved pipeline output into stage subfolders
# (Catalog/ Calls/ Haps/ Scans/); projects predating it have RefAlt.*, R.haps.*, and
# <scan>/ flat in process/<project>/. This reports each stage in BOTH places, so a
# half-migrated project is visible rather than looking empty.
#
# "layout: MIXED" is the one to act on: a project with RefAlt both at top level and
# in Calls/ must NOT be fed to reorganize_project.sh as-is — it does a bare
# `mv RefAlt.*.txt Calls/`, so top-level (often stale symlinks) overwrite the real
# merged tables in Calls/. Remove the redundant top-level copies first.
set -uo pipefail
pat="${1:?usage: process_fingerprint.sh <prefix>   e.g. AGE_SY, ZINC}"
OUT="logs/fingerprint.${pat}.txt"; mkdir -p logs

# count files matching a glob in one dir, no recursion; 0 if the dir is absent
count_in () {  # count_in <dir> <name-pattern> [-type l|f]
  find "$1" -maxdepth 1 -name "$2" "${@:3}" 2>/dev/null | wc -l | tr -d ' '
}
# space-separated basenames of subdirs, excluding the four stage dirs
dirs_in () {
  find "$1" -maxdepth 1 -mindepth 1 -type d 2>/dev/null \
    | sed 's|.*/||' | grep -vxE 'Catalog|Calls|Haps|Scans' | sort | tr '\n' ' '
}

{
shopt -s nullglob
for d in process/${pat}*/; do
  d=${d%/}
  echo "################ $d ################"
  echo "  newest mtime : $(find "$d" -maxdepth 4 -type f -printf '%TY-%Tm-%Td %TH:%TM\n' 2>/dev/null | sort | tail -1)"

  # ── RefAlt: top-level (pre-#28) vs Calls/ (current) ───────────────────────
  ra_top_l=$(count_in "$d" 'RefAlt.chr*.txt' -type l)
  ra_top_f=$(count_in "$d" 'RefAlt.chr*.txt' -type f)
  ra_top=$((ra_top_l + ra_top_f))
  ra_calls=$(count_in "$d/Calls" 'RefAlt.chr*.txt')

  # ── haplotypes / scans in both places ────────────────────────────────────
  haps_top=$(count_in "$d" 'R.haps.*.rds')
  haps_new=$(count_in "$d/Haps" 'R.haps.*.rds')
  scans_top=$(dirs_in "$d")
  scans_new=$(dirs_in "$d/Scans")

  # ── layout verdict ───────────────────────────────────────────────────────
  staged=0; flat=0
  [ "$ra_calls" -gt 0 ] || [ "$haps_new" -gt 0 ] || [ -n "$scans_new" ] && staged=1
  [ "$ra_top"   -gt 0 ] || [ "$haps_top" -gt 0 ] || [ -n "$scans_top" ] && flat=1
  if   [ "$staged" = 1 ] && [ "$flat" = 1 ]; then layout="MIXED — migrate with care (see header)"
  elif [ "$staged" = 1 ];                    then layout="staged (Calls/Haps/Scans)"
  elif [ "$flat"   = 1 ];                    then layout="flat (pre-#28) — needs reorganize_project.sh"
  else                                            layout="empty / nothing recognized"
  fi
  echo "  layout       : $layout"
  if [ "$ra_top" -gt 0 ] && [ "$ra_calls" -gt 0 ]; then
    echo "  !! CLOBBER RISK: RefAlt exists at top level ($ra_top) AND in Calls/ ($ra_calls)."
    echo "  !! reorganize_project.sh would overwrite Calls/ with the top-level copies."
  fi

  echo "  Catalog/     : $([ -d "$d/Catalog" ] && echo "yes ($(count_in "$d/Catalog" '*') files)" || echo no)"
  echo "  calls.bcf    : top=$(count_in "$d" 'calls.*.bcf')  Calls/=$(count_in "$d/Calls" 'calls.*.bcf')"
  echo "  counts/      : $([ -d "$d/Calls/counts" ] && echo "Calls/counts ($(count_in "$d/Calls/counts" '*.tsv.gz'))" \
                          || { [ -d "$d/counts" ] && echo "counts/ (top, $(count_in "$d/counts" '*.tsv.gz'))" || echo none; })"
  echo "  RefAlt       : top=$ra_top (real=$ra_top_f symlink=$ra_top_l)  Calls/=$ra_calls"
  echo "  RefAlt/founder (fromsites): $(find "$d" "$d/Calls" -maxdepth 1 -name 'RefAlt.*.chr*.txt' 2>/dev/null | grep -cvE '/RefAlt\.chr')"
  echo "  R.haps       : top=$haps_top  Haps/=$haps_new"
  echo "  scan dirs    : top=[${scans_top}] Scans/=[${scans_new}]"

  # samples from a whole-sample RefAlt, wherever it lives (resolves symlinks)
  ra=$(ls "$d"/Calls/RefAlt.chrX.txt "$d"/RefAlt.chrX.txt 2>/dev/null | head -1)
  if [ -n "$ra" ]; then
    samps=$(head -1 "$ra" 2>/dev/null | tr ' \t' '\n' | sed -n 's/^REF_//p')
    n=$(printf '%s\n' "$samps" | grep -c .)
    echo "  samples ($n): $(printf '%s\n' "$samps" | head -8 | tr '\n' ' ')..."
  fi
  echo
done
} | tee "$OUT"
echo "wrote $OUT  ->  bash scripts/cluster_sync.sh"
