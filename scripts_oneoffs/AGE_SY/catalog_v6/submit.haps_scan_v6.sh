#!/bin/bash
# submit.haps_scan_v6.sh — Steps 3+4 chained, one command:
#   bash scripts_oneoffs/AGE_SY/catalog_v6/submit.haps_scan_v6.sh
#
# Step 3: call reps 7-12 against the v6 catalog -> re-merge v6/Calls/RefAlt.* to R1-R12.
# Step 4: size75k haplotypes on the v6 RefAlt -> the 4 scans (smooth 100) -> the 4-scan
#         figure -> figures/AGE_SY_v6_size75k_4scan.png. Reuses the round3 run_scans.sh +
#         plot_4scan.sh (both env-parameterized) unchanged.
# Everything chains by SLURM dependency. Compare the result to the July-15
# figures/AGE_SY_v3_size75k_4scan.png (both copied into logs/figures/ at the end).
set -euo pipefail
CATDIR=process/AGE_SY_v6/Catalog
V6=process/AGE_SY_v6
SRC=process/AGE_SY_v6/Calls                 # where v6 RefAlt.<chr>.txt live
NEW=process/AGE_SY_v6_size75k               # size75k haps+scan output (parallels AGE_SY_v3_size75k)
BAMLIST=helpfiles/AGE_SY/bam_list.v4_R7-12.txt
PAR_SRC=helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R
PAR_NEW=helpfiles/AGE_SY/AGE_SY_haplotype_parameters_size75k.R
SIZE=75000; SMOOTH=100
CHRS="chrX chr2L chr2R chr3L chr3R"
R3=scripts_oneoffs/AGE_SY/round3_v3_R12
for f in "$BAMLIST" "$PAR_SRC" "$R3/run_scans.sh" "$R3/plot_4scan.sh" "$CATDIR/catalog.tsv.gz"; do
  [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

# --- Step 3: call reps 7-12 (adds to v6/Calls/counts, re-merges RefAlt to R1-R12) ---
JID_CALL=$(bash pipeline/scripts/call_samples.sh --catalog "$CATDIR" --bamlist "$BAMLIST" --dir "$V6" | tail -1)
echo "Step 3 call R7-12 (merge job): $JID_CALL"
[[ "$JID_CALL" =~ ^[0-9]+$ ]] || { echo "ERROR: no job id from call_samples (got '$JID_CALL')" >&2; exit 1; }

# --- Step 4a: size75k parameter file (one line changed from the real one) ---
sed 's/^size *=.*/size = '"$SIZE"'/' "$PAR_SRC" > "$PAR_NEW"
echo "parfile: $PAR_NEW"; grep -E '^(step|size) *=' "$PAR_NEW" | sed 's/^/  /'

# --- Step 4b: haps dir; symlink the v6 RefAlt in (content filled by the R7-12 merge) ---
mkdir -p "$NEW"
for c in $CHRS; do ln -sf "$(cd "$SRC" && pwd)/RefAlt.$c.txt" "$NEW/RefAlt.$c.txt"; done

# --- Step 4c: haplotypes (size75k) AFTER the R7-12 merge ---
JID_HAPS=$(sbatch --parsable --dependency=afterok:"$JID_CALL" \
    pipeline/scripts/REFALT2haps.sh --parfile "$PAR_NEW" --dir "$NEW")
echo "haps size75k (after R7-12): $JID_HAPS -> $NEW"

# --- Step 4d: the 4 scans (smooth 100) after haps (reuse round3 run_scans.sh) ---
scan_out=$(AFTER="$JID_HAPS" SMOOTH="$SMOOTH" DIR="$NEW" bash "$R3/run_scans.sh")
printf '%s\n' "$scan_out"
CONCAT=$(printf '%s\n' "$scan_out" | sed -n 's/^CONCAT_JIDS=//p')
[[ -n "$CONCAT" ]] || { echo "ERROR: no scan concat job ids" >&2; exit 1; }

# --- Step 4e: the figure after the scans (reuse round3 plot_4scan.sh) ---
mkdir -p figures
JID_FIG=$(sbatch --parsable --dependency=afterok:"$CONCAT" --export=ALL,SCAN_DIR="$NEW" "$R3/plot_4scan.sh")
echo "figure: $JID_FIG -> figures/AGE_SY_v6_size75k_4scan.png"

# --- final: copy v6 + the July-15 v3 figure into logs/figures/ so both sync to Claude ---
sbatch --parsable --dependency=afterok:"$JID_FIG" -A tdlong_lab -p standard --time=5:00 \
  -o logs/AGE_SY/copy_scanfig.out --wrap="mkdir -p logs/figures; \
     cp -f figures/AGE_SY_v6_size75k_4scan.png figures/AGE_SY_v3_size75k_4scan.png logs/figures/ 2>/dev/null; echo copied" >/dev/null
echo
echo "chain: call R7-12 ($JID_CALL) -> haps ($JID_HAPS) -> 4 scans -> figure ($JID_FIG) -> copy to logs/figures/"
squeue -u "$USER"
