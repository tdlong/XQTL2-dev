#!/bin/bash
###############################################################################
# 03_update_helpfiles.sh — AGE_SY round 3 (v3): regenerate bam list + params
#
# Run AFTER 01_align.sh finishes. Rebuilds, from what is actually on disk:
#   - helpfiles/AGE_SY/AGE_SY.bams               (sample BAMs >1 GB + founders)
#   - helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R  (names_in_bam line)
#
# Both are written in the project's canonical order:
#   for rep 1..12: for sex F,M: for trt AgeSY10,AgeSY20,Con
# and the two files are guaranteed consistent (same samples, same order).
#
# BAMs <1 GB (or missing) are excluded automatically — the pipeline drops
# low-coverage libraries, so they must not appear in either file.
#
# The founder block at the bottom of AGE_SY.bams is preserved verbatim.
#
# NOT handled here (needs real data — the only manual step): the four
# design files. Add R12 fly counts to summary_info_v1.xlsx, then run
#   python3 scripts_oneoffs/AGE_SY/common/make_AGE_SY_design_files.py
#
# Run from: XQTL2-dev repo root, on the cluster.
###############################################################################
set -euo pipefail

BAMDIR=data/bam/AGE_SY
BAMLIST=helpfiles/AGE_SY/AGE_SY.bams
PARFILE=helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R
REPS=$(seq 1 12)
SEXES="F M"
TRTS="AgeSY10 AgeSY20 Con"

[ -d "$BAMDIR" ]  || { echo "ERROR: $BAMDIR not found (run on the cluster)" >&2; exit 1; }
[ -f "$BAMLIST" ] || { echo "ERROR: $BAMLIST not found" >&2; exit 1; }
[ -f "$PARFILE" ] || { echo "ERROR: $PARFILE not found" >&2; exit 1; }

# --- collect samples present and >1 GB, in canonical order ---
samples=()
for rep in $REPS; do for sex in $SEXES; do for trt in $TRTS; do
    s="${trt}_R${rep}_${sex}"; f="${BAMDIR}/${s}.bam"
    if [ -f "$f" ]; then
        if [ -n "$(find "$f" -size +1G 2>/dev/null)" ]; then
            samples+=("$s")
        else
            echo "  skip (<1 GB): $s"
        fi
    fi
done; done; done

[ "${#samples[@]}" -gt 0 ] || { echo "ERROR: no sample BAMs >1 GB found" >&2; exit 1; }

# sanity: warn about any *.bam >1 GB in the dir that the grid did not match
grid_n=${#samples[@]}
disk_n=$(find "$BAMDIR" -maxdepth 1 -name '*.bam' -size +1G | wc -l | tr -d ' ')
if [ "$grid_n" != "$disk_n" ]; then
    echo "WARNING: ${disk_n} BAMs >1 GB on disk but ${grid_n} matched the naming grid" >&2
    echo "         (a BAM may not fit {AgeSY10,AgeSY20,Con}_R{1..12}_{F,M}.bam)" >&2
fi

# --- rewrite AGE_SY.bams: samples (relative paths) + preserved founder block ---
founders=$(grep -v "^${BAMDIR}/" "$BAMLIST")        # the non-sample lines = founders
{
    for s in "${samples[@]}"; do echo "${BAMDIR}/${s}.bam"; done
    echo "$founders"
} > "$BAMLIST"

# --- rewrite the names_in_bam line in the params file (leave everything else) ---
names=$(printf '"%s",' "${samples[@]}"); names="${names%,}"   # "a","b",...
grep -v '^names_in_bam' "$PARFILE" > "${PARFILE}.tmp"
echo "names_in_bam=c(${names})" >> "${PARFILE}.tmp"
mv "${PARFILE}.tmp" "$PARFILE"

echo ""
echo "Wrote ${#samples[@]} samples to:"
echo "  $BAMLIST  (+ $(echo "$founders" | grep -c . ) founders)"
echo "  $PARFILE  (names_in_bam)"
echo ""
echo "R12 samples included:"
printf '  %s\n' "${samples[@]}" | grep _R12_ || echo "  (none — check R12 BAM sizes)"
echo ""
echo "Next: update the 4 design files (add R12 to summary_info_v1.xlsx, then"
echo "  python3 scripts_oneoffs/AGE_SY/common/make_AGE_SY_design_files.py),"
echo "  commit + push, then run 02_refalt_haps_scan.sh"
