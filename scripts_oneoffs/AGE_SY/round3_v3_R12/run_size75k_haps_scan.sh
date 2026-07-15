#!/bin/bash
###############################################################################
# run_size75k_haps_scan.sh — re-run the haplotype step at HALF the window
# (size=75000 instead of 150000) into a NEW folder, then the four scans at
# smooth=100, so we can compare against the 150 kb version.
#
#   process/AGE_SY_v3        <- untouched, 150 kb haplotype window
#   process/AGE_SY_v3_size75k <- new, 75 kb haplotype window
#
# 'size' is a config switch (in the haplotype parameter file), not a code edit:
# this derives a size=75000 parfile from the real one by changing one line.
# Haplotype estimation is re-done (the RefAlt SNP tables are reused via symlink),
# then smooth(100) -> scan -> concat, all into the new folder.
#
# Run once from repo root on the cluster:
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/run_size75k_haps_scan.sh
###############################################################################
set -euo pipefail
HERE=scripts_oneoffs/AGE_SY/round3_v3_R12
SRC=process/AGE_SY_v3
NEW=process/AGE_SY_v3_size75k
PAR_SRC=helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R
PAR_NEW=helpfiles/AGE_SY/AGE_SY_haplotype_parameters_size75k.R
SIZE=75000
SMOOTH=100
CHRS="chrX chr2L chr2R chr3L chr3R"

# guards: the RefAlt SNP tables must exist (we reuse them, not recompute)
for c in $CHRS; do
  [ -f "$SRC/RefAlt.$c.txt" ] || { echo "ERROR: missing $SRC/RefAlt.$c.txt" >&2; exit 1; }
done

# 1. new folder; symlink the RefAlt tables in (REFALT2haps reads them from --dir)
mkdir -p "$NEW"
for c in $CHRS; do ln -sf "$(cd "$SRC" && pwd)/RefAlt.$c.txt" "$NEW/RefAlt.$c.txt"; done

# 2. size=75000 parameter file, derived from the real one (one line changed)
sed 's/^size *=.*/size = '"$SIZE"'/' "$PAR_SRC" > "$PAR_NEW"
echo "parameter file: $PAR_NEW"
grep -E '^(step|size) *=' "$PAR_NEW" | sed 's/^/  /'

# 3. haplotypes at size=75000 into the new folder (direct call = 1-day budget)
JID_HAPS=$(sbatch --parsable pipeline/scripts/REFALT2haps.sh \
    --parfile "$PAR_NEW" --dir "$NEW")
echo "haplotypes (size=$SIZE) submitted: $JID_HAPS  -> $NEW"

# 4. the four scans at smooth=100 into the new folder, after haps completes
scan_out=$(AFTER="$JID_HAPS" SMOOTH="$SMOOTH" DIR="$NEW" bash "$HERE/run_scans.sh")
printf '%s\n' "$scan_out"
CONCAT=$(printf '%s\n' "$scan_out" | sed -n 's/^CONCAT_JIDS=//p')
[ -n "$CONCAT" ] || { echo "ERROR: could not capture scan concat job ids" >&2; exit 1; }

# 5. the figure, after all four scans finish (no manual step)
mkdir -p figures
JID_FIG=$(sbatch --parsable --dependency=afterok:"$CONCAT" \
    --export=ALL,SCAN_DIR="$NEW" "$HERE/plot_4scan.sh")
echo "figure submitted (after the 4 scans): $JID_FIG"

cat <<EOF

================================================================
Queued as ONE chain — nothing else to run (150 kb version in $SRC untouched):
  haps size=$SIZE ($JID_HAPS) -> smooth ${SMOOTH}kb -> scan -> concat -> figure ($JID_FIG)
Monitor: squeue -u \$USER
Final figure appears at: figures/$(basename "$NEW")_4scan.png
  (compare to figures/AGE_SY_v3_4scan.png)
================================================================
EOF
