#!/bin/bash
###############################################################################
# inventory.sh — what malathion data actually exists on HPC3, returned via git.
#
# Read-only, runs in well under a minute, fine on the login node. Answers the
# one question the manifests cannot: which pools are in the RefAlt tables
# (rep 1 only, or all 8 of the 2022 paper), and what the Nov 2023 reanalysis
# directory behind process/MALATHION_v1 holds. That directory is outside the
# repo, so it is reached by resolving the existing symlink, never by path.
#
#   bash scripts_oneoffs/malathion/inventory.sh
#   bash scripts/cluster_sync.sh
#
# Writes logs/malathion/inventory.txt (and the fingerprint files).
###############################################################################
set -uo pipefail
[ -d scripts_oneoffs/malathion ] || { echo "run from the XQTL2-dev root" >&2; exit 1; }
mkdir -p logs/malathion
OUT=logs/malathion/inventory.txt

{
echo "## $(date)  $(hostname)"
echo
echo "## process fingerprints"
bash scripts/process_fingerprint.sh MALATHION | grep -v '^wrote'
bash scripts/process_fingerprint.sh malathion | grep -v '^wrote'
bash scripts/process_fingerprint.sh XQTL_talk | grep -v '^wrote'

echo "## RefAlt sample columns (the whole header, one dir per block)"
for d in process/MALATHION_v1 process/XQTL_talk/MALATHION process/malathion_test; do
  f=$d/RefAlt.chrX.txt
  [ -e "$f" ] || { echo "-- $f: missing"; continue; }
  echo "-- $f"
  echo "   resolves to: $(readlink -f "$f")"
  echo "   size: $(du -h "$(readlink -f "$f")" | cut -f1)   rows: $(wc -l < "$f")"
  echo "   columns: $(head -1 "$f" | tr ' \t' '\n' | grep -v '^$' | tr '\n' ' ')"
done
echo

echo "## samples inside the existing haplotype fits"
module load R/4.2.2 2>/dev/null
for d in process/MALATHION_v1 process/malathion_test; do
  f=$d/R.haps.chrX.out.rds
  [ -e "$f" ] || { echo "-- $f: missing"; continue; }
  echo "-- $f"
  Rscript -e 'x <- readRDS(commandArgs(TRUE)[1]); cat("   windows:", nrow(x), " samples:", paste(x$sample[[1]], collapse=" "), "\n")' "$f"
done
echo

echo "## the Nov 2023 reanalysis directory (resolved from the MALATHION_v1 symlink)"
base=$(dirname "$(dirname "$(readlink -f process/MALATHION_v1/RefAlt.chrX.txt)")")
echo "-- base: $base"
ls -la "$base" 2>/dev/null
echo "-- BAMs:"
ls -la "$base"/data/bam 2>/dev/null || echo "   (no data/bam)"
echo "-- helperfiles:"
for f in "$base"/helperfiles/*; do
  [ -e "$f" ] || continue
  echo "== $f"
  head -20 "$f"
done
echo

echo "## BAMs in this repo"
ls -la data/bam/malathion* 2>/dev/null || echo "   (none)"
} 2>&1 | tee "$OUT"

echo
echo "wrote $OUT  ->  bash scripts/cluster_sync.sh"
