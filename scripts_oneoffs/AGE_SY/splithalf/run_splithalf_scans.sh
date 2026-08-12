#!/bin/bash
###############################################################################
# run_splithalf_scans.sh — the eight split-half scans, into their own folder.
#
# This is NOT the main pipeline. It writes only to process/AGE_SY_splithalf and
# never inside process/AGE_SY: the haplotypes are reused read-only through a
# symlink, so the live scans and figures cannot be touched or overwritten.
#
#   process/AGE_SY            <- live, untouched
#   process/AGE_SY_splithalf  <- everything this produces
#     Haps -> ../AGE_SY/Haps    (symlink; run_scan.sh reads <dir>/Haps)
#     Scans/AGE_SY10_F_odd/ ... (8 scans, written under <dir>/Scans/<scan>)
#
# smooth=100 kb to match the live run (haplotypes are size=75000).
#
# Run once from repo root on the cluster:
#   bash scripts_oneoffs/AGE_SY/splithalf/run_splithalf_scans.sh
###############################################################################
set -euo pipefail

LIVE=process/AGE_SY
NEW=process/AGE_SY_splithalf
DESIGN_DIR=helpfiles/AGE_SY/splithalf
SMOOTH=100

# ── guards ───────────────────────────────────────────────────────────────────
[ -d "$LIVE/Haps" ] || { echo "ERROR: missing $LIVE/Haps" >&2; exit 1; }
for scan in AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M; do
  for half in odd even; do
    d="$DESIGN_DIR/$scan.$half.txt"
    [ -f "$d" ] || { echo "ERROR: missing $d — run make_splithalf_designs.R first" >&2; exit 1; }
  done
done

# ── output folder, haplotypes symlinked in read-only ─────────────────────────
mkdir -p "$NEW"
[ -e "$NEW/Haps" ] || ln -s "$(cd "$LIVE" && pwd)/Haps" "$NEW/Haps"
echo "output:     $NEW"
echo "haplotypes: $(readlink "$NEW/Haps")  (reused, not recomputed)"

# ── the eight scans ──────────────────────────────────────────────────────────
JIDS=""
for scan in AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M; do
  for half in odd even; do
    name="${scan}_${half}"
    out=$(bash pipeline/scripts/run_scan.sh \
            --smooth "$SMOOTH" --dir "$NEW" --scan "$name" \
            --design "$DESIGN_DIR/$scan.$half.txt")
    jid=$(printf '%s\n' "$out" | awk '/^done:/{print $2}')
    [ -n "$jid" ] || { echo "ERROR: no concat job id for $name" >&2; exit 1; }
    JIDS="${JIDS:+$JIDS:}$jid"
    echo "submitted:  $name (smooth ${SMOOTH}kb, concat job $jid)"
  done
done

# ── gather the eight Cutl_H2 columns into one long file, after all finish ────
mkdir -p logs/AGE_SY
JID_GATHER=$(sbatch --parsable --dependency=afterok:"$JIDS" \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=6G --time=30:00 \
    -o logs/AGE_SY/splithalf_gather.%j.out \
    --wrap="module load R/4.2.2; Rscript scripts_oneoffs/AGE_SY/splithalf/gather_splithalf_H2.R")
echo "gather:     $JID_GATHER (after all 8 scans)"

cat <<EOF

All 8 scans submitted -> $NEW/Scans/<scan>_<half>/
Gather job $JID_GATHER writes $NEW/AGE_SY_splithalf_H2.txt.gz
EOF
