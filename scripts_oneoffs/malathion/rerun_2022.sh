#!/bin/bash
###############################################################################
# rerun_2022.sh — the full 2022 malathion experiment through the current
# pipeline, from RefAlt, into its own project folder.
#
# Why: Stuart's rerun of the published data gives chr2R ~20 and chr3L ~5 where
# the 2022 paper had both near 5. The suspected cause is the test, not the data
# (see temp_malathion/). This reproduces his run here and adds what diagnoses it:
#
#   <prefix>_all   every pool in the design, all REP   (Stuart's run, if 8 pools)
#   <prefix>_R<r>  one scan per REP label in the design: a single control/selected
#                  pair. The Wald test weights pairs by pool size, so these show
#                  what each pair contributes on its own.
#   finish_2022.sh the paper's ANOVA on the _all haplotypes, Fig.4-style shift
#                  tables at the two peaks, a peak summary -> logs/malathion/
#
# The design decides how many pairs there are. design_2022.txt is the paper's
# 8 pools (REP 1-4); design_rep1.txt is the 4 rep-1 pools that exist on HPC3
# (REP 1-2). For each REP label r the file <design>_R<r>.txt must exist.
#
# Nothing is written inside process/MALATHION_v1 or any other existing folder.
# RefAlt is reached through relative symlinks; haplotypes are recomputed with
# --parfile unless --reuse-haps points at fits of the same pools.
#
# Before submitting it checks the RefAlt header for every pool in the design
# and every founder, and stops if one is missing.
#
# HPC3 holds rep 1 only (logs/malathion/inventory.txt, 2026-11-24), so here:
#   bash scripts_oneoffs/malathion/rerun_2022.sh \
#       --design helpfiles/malathion/design_rep1.txt \
#       --parfile helpfiles/malathion/hap_params.R \
#       --reuse-haps process/MALATHION_v1 \
#       --dir process/malathion_rep1 --prefix malathion_rep1
# With all 8 pools (Stuart's data) the defaults apply and haplotypes are refit.
# Options:
#   --refalt <dir>      folder holding RefAlt.chr*.txt   (default process/MALATHION_v1)
#   --dir <dir>         project folder to create         (default process/malathion_2022)
#   --prefix <name>     scan name prefix                 (default malathion_2022)
#   --design <file>     design, any number of REP        (default helpfiles/malathion/design_2022.txt)
#   --parfile <file>    haplotype parameters             (default helpfiles/malathion/hap_params_2022.R)
#   --smooth <kb>       smoothing half-window            (default 250, the pipeline default)
#   --reuse-haps <dir>  folder holding R.haps.chr*.out.rds fit from the SAME pools;
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
# REP labels in the design (write.table layout: rowname bam TRT REP ...)
REPS=$(awk 'NR>1 {gsub(/"/,"",$4); print $4}' "$DESIGN" | sort -nu)
[ -n "$REPS" ] || { echo "ERROR: no REP column values in $DESIGN" >&2; exit 1; }
for r in $REPS; do
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

# ── the scans: all pools, then one per pair ──────────────────────────────────
# Single-sex pairs get --sex so chrX is credited the right number of X's. The
# all-pools scan mixes sexes across pairs; chrX there is not interpretable and
# is left at the default. Autosomes are unaffected either way.
submit_scan () {  # submit_scan <name> <design> <sex>
  local out jid
  out=$(bash pipeline/scripts/run_scan.sh --design "$2" --dir "$DIR" --scan "$1" \
          --smooth "$SMOOTH" --sex "$3" ${AFTER[@]+"${AFTER[@]}"})
  jid=$(printf '%s\n' "$out" | awk '/^done:/{print $2}')
  [ -n "$jid" ] || { echo "ERROR: no concat job id for $1" >&2; printf '%s\n' "$out" >&2; exit 1; }
  JIDS="${JIDS:+$JIDS:}$jid"
  echo "scan:       $1  ($2, --sex $3, concat job $jid)"
}
JIDS=""
ALSO=""
submit_scan "${PREFIX}_all" "$DESIGN" mixed
for r in $REPS; do
  dr="$DESIGN_DIR/${DESIGN_BASE}_R$r.txt"
  # sex of the pair from its sample names (c_F_1 / s_M_2 ...); mixed if unclear
  sex=mixed
  if   grep -q '_F_' "$dr" && ! grep -q '_M_' "$dr"; then sex=F
  elif grep -q '_M_' "$dr" && ! grep -q '_F_' "$dr"; then sex=M; fi
  submit_scan "${PREFIX}_R$r" "$dr" "$sex"
  ALSO="${ALSO:+$ALSO,}${PREFIX}_R$r"
done

# ── diagnostics + return through git, after every scan finishes ─────────────
JID_FIN=$(sbatch --parsable --dependency=afterok:"$JIDS" \
    scripts_oneoffs/malathion/finish_2022.sh "$DIR" "$PREFIX" "$ALSO")
echo "finish:     job $JID_FIN  (paper ANOVA, peak tables, figures -> logs/malathion/, then git push)"
echo
echo "When it is done, on the laptop: git pull, then read logs/malathion/."
