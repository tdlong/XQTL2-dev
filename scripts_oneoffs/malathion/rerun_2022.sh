#!/bin/bash
###############################################################################
# rerun_2022.sh — the full 2022 malathion experiment through the current
# pipeline, from RefAlt, into its own project folder.
#
# Why: Stuart's rerun of the published data gives chr2R ~20 and chr3L ~5 where
# the 2022 paper had both near 5. The suspected cause is the test, not the data
# (see temp_malathion/). This reproduces his run here and adds what diagnoses it:
#
#   <prefix>_all   all 8 pools, REP 1-4           (Stuart's run)
#   <prefix>_R1    F rep 1 only   (43 vs 43)       single-pair scans: the Wald
#   <prefix>_R2    M rep 1 only   (65 vs 65)       test weights pairs by pool
#   <prefix>_R3    F rep 2 only  (264 vs 693)      size, so these show what
#   <prefix>_R4    M rep 2 only  (264 vs 403)      each pair contributes
#   finish_2022.sh the paper's ANOVA on the _all haplotypes, Fig.4-style shift
#                  tables at the two peaks, a peak summary -> logs/malathion/
#
# Nothing is written inside process/MALATHION_v1 or any other existing folder.
# RefAlt is reached through relative symlinks; haplotypes are recomputed with
# helpfiles/malathion/hap_params_2022.R (8 pools) unless --reuse-haps is given.
#
# Before submitting it checks the RefAlt header for every pool in the design
# and every founder, and stops if one is missing -- the rep-2 sample names in
# the design are the rep-1 pattern and have not been seen in a header yet.
#
# Run once, from the repo root on HPC3:
#   bash scripts_oneoffs/malathion/rerun_2022.sh
# Options:
#   --refalt <dir>      folder holding RefAlt.chr*.txt   (default process/MALATHION_v1)
#   --dir <dir>         project folder to create         (default process/malathion_2022)
#   --prefix <name>     scan name prefix                 (default malathion_2022)
#   --design <file>     8-pool design                    (default helpfiles/malathion/design_2022.txt)
#   --parfile <file>    haplotype parameters             (default helpfiles/malathion/hap_params_2022.R)
#   --smooth <kb>       smoothing half-window            (default 250, the pipeline default)
#   --reuse-haps <dir>  folder holding R.haps.chr*.out.rds fit from the SAME 8 pools;
#                       symlinked instead of recomputed (skips a day of highmem)
###############################################################################
set -euo pipefail
[ -d scripts_oneoffs/malathion ] || { echo "ERROR: run from the XQTL2-dev root" >&2; exit 1; }

REFALT=process/MALATHION_v1
DIR=process/malathion_2022
PREFIX=malathion_2022
DESIGN=helpfiles/malathion/design_2022.txt
PARFILE=helpfiles/malathion/hap_params_2022.R
SMOOTH=250
REUSE_HAPS=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --refalt)     REFALT="$2";     shift 2 ;;
    --dir)        DIR="$2";        shift 2 ;;
    --prefix)     PREFIX="$2";     shift 2 ;;
    --design)     DESIGN="$2";     shift 2 ;;
    --parfile)    PARFILE="$2";    shift 2 ;;
    --smooth)     SMOOTH="$2";     shift 2 ;;
    --reuse-haps) REUSE_HAPS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done
DESIGN_DIR=$(dirname "$DESIGN")
DESIGN_BASE=$(basename "$DESIGN" .txt)
chrs=(chrX chr2L chr2R chr3L chr3R)

# ── guards ───────────────────────────────────────────────────────────────────
for chr in "${chrs[@]}"; do
  [ -e "$REFALT/RefAlt.$chr.txt" ] || { echo "ERROR: missing $REFALT/RefAlt.$chr.txt" >&2; exit 1; }
done
[ -f "$DESIGN" ]  || { echo "ERROR: missing $DESIGN" >&2; exit 1; }
[ -f "$PARFILE" ] || { echo "ERROR: missing $PARFILE" >&2; exit 1; }
for r in 1 2 3 4; do
  [ -f "$DESIGN_DIR/${DESIGN_BASE}_R$r.txt" ] || { echo "ERROR: missing $DESIGN_DIR/${DESIGN_BASE}_R$r.txt" >&2; exit 1; }
done
if [ -e "$DIR/Scans/${PREFIX}_all" ]; then
  echo "ERROR: $DIR/Scans/${PREFIX}_all exists; pick another --prefix or --dir rather than overwrite" >&2; exit 1
fi

# Every pool in the design, and every founder, must be a column of RefAlt.
header=$(head -1 "$REFALT/RefAlt.chrX.txt" | tr ' \t' '\n' | grep -v '^$')
pools=$(awk 'NR>1 {gsub(/"/,"",$2); print $2}' "$DESIGN")
founders=$(Rscript -e 'source(commandArgs(TRUE)[1]); cat(founders)' "$PARFILE")
missing=""
for s in $pools $founders; do
  printf '%s\n' "$header" | grep -qx "REF_$s" || missing="$missing $s"
done
if [ -n "$missing" ]; then
  echo "ERROR: not in $REFALT/RefAlt.chrX.txt header:$missing" >&2
  echo "       header samples: $(printf '%s\n' "$header" | sed -n 's/^REF_//p' | tr '\n' ' ')" >&2
  echo "       If rep 2 is absent here, the 8-pool run needs a RefAlt that has it (--refalt)." >&2
  exit 1
fi
echo "RefAlt:     $REFALT  (all $(echo $pools | wc -w) pools and $(echo $founders | wc -w) founders present)"

# ── project folder: RefAlt linked in read-only, stage layout ─────────────────
mkdir -p "$DIR/Calls" "$DIR/Haps" "$DIR/Scans" logs/malathion
for chr in "${chrs[@]}"; do
  dst="$DIR/Calls/RefAlt.$chr.txt"
  [ -e "$dst" ] || ln -s "$(realpath -s --relative-to="$DIR/Calls" "$REFALT/RefAlt.$chr.txt")" "$dst"
done
echo "project:    $DIR  (Calls/RefAlt -> $(readlink "$DIR/Calls/RefAlt.chrX.txt"))"

# ── haplotypes ───────────────────────────────────────────────────────────────
if [ -n "$REUSE_HAPS" ]; then
  for chr in "${chrs[@]}"; do
    src="$REUSE_HAPS/R.haps.$chr.out.rds"
    [ -e "$src" ] || { echo "ERROR: missing $src" >&2; exit 1; }
    dst="$DIR/Haps/R.haps.$chr.out.rds"
    [ -e "$dst" ] || ln -s "$(realpath -s --relative-to="$DIR/Haps" "$src")" "$dst"
  done
  AFTER=()
  echo "haplotypes: reused from $REUSE_HAPS (not recomputed)"
else
  JID_HAPS=$(bash pipeline/scripts/run_haps.sh --parfile "$PARFILE" --dir "$DIR")
  AFTER=(--after "$JID_HAPS")
  echo "haplotypes: job $JID_HAPS  (REFALT2haps, 5 chromosomes, highmem)"
fi

# ── the five scans ───────────────────────────────────────────────────────────
# Single-sex pairs get --sex so chrX is credited the right number of X's. The
# all-pools scan mixes sexes across pairs; chrX there is not interpretable and
# is left at the default. Autosomes are unaffected either way.
submit_scan () {  # submit_scan <name> <design> <sex>
  local out jid
  out=$(bash pipeline/scripts/run_scan.sh --design "$2" --dir "$DIR" --scan "$1" \
          --smooth "$SMOOTH" --sex "$3" "${AFTER[@]}")
  jid=$(printf '%s\n' "$out" | awk '/^done:/{print $2}')
  [ -n "$jid" ] || { echo "ERROR: no concat job id for $1" >&2; printf '%s\n' "$out" >&2; exit 1; }
  JIDS="${JIDS:+$JIDS:}$jid"
  echo "scan:       $1  ($2, --sex $3, concat job $jid)"
}
JIDS=""
submit_scan "${PREFIX}_all" "$DESIGN" mixed
submit_scan "${PREFIX}_R1"  "$DESIGN_DIR/${DESIGN_BASE}_R1.txt" F
submit_scan "${PREFIX}_R2"  "$DESIGN_DIR/${DESIGN_BASE}_R2.txt" M
submit_scan "${PREFIX}_R3"  "$DESIGN_DIR/${DESIGN_BASE}_R3.txt" F
submit_scan "${PREFIX}_R4"  "$DESIGN_DIR/${DESIGN_BASE}_R4.txt" M

# ── diagnostics + return through git, after all five finish ─────────────────
JID_FIN=$(sbatch --parsable --dependency=afterok:"$JIDS" \
    scripts_oneoffs/malathion/finish_2022.sh "$DIR" "$PREFIX")
echo "finish:     job $JID_FIN  (paper ANOVA, peak tables, figures -> logs/malathion/, then git push)"
echo
echo "When it is done, on the laptop: git pull, then read logs/malathion/."
