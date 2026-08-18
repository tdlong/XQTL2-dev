#!/bin/bash
# run_all.sh — put the Aug 2024 aging pilot through the current pipeline, and
# re-scan AGE_SY on its first 6 replicates so the two can be compared.
#
# Run from the repo root on HPC3, after `git pull`:
#   bash scripts_oneoffs/AGE_2024/run_all.sh
#
# Submits everything in one go and prints the job IDs. Nothing here runs the
# analysis itself -- it chains the pipeline's own scripts.
#
# WHAT RUNS, AND FROM WHAT
#
#   1. copy the SNP list from AGE_SY        process/AGE_SY/Catalog -> process/AGE_2024
#      Same 8 founders, so the same SNP list applies, and copying rather than
#      rebuilding means both experiments are counted at identical sites. The 8
#      founder count files come across too: same BAMs against the same catalog,
#      so they are identical by construction and need not be recomputed.
#
#   2. call SNPs per sample                 12 BAMs in data/bam/AGE_2024
#      One job per sample, then a merge into RefAlt.<chr>.txt that picks up the
#      copied founder counts alongside the new sample counts.
#
#   3. call haplotypes                      75 kb window, 5 kb step
#      Same parameters as AGE_SY's size75k run, so the scans are comparable.
#
#   4. scan                                 smoothed at 100 kb
#      Design is helpfiles/AGE_Aug13_24/Ageing_Aug13.txt, already in the repo:
#      6 cages, one control and one selected pool each, females on lab food.
#
#   5. re-scan AGE_SY on reps 1-6           4 scans, into process/AGE_SY
#      Haplotypes already exist there, so this is scan-only and cheap. Runs
#      immediately, in parallel with steps 2-4.
#
# Steps 2-4 are chained on job IDs. Step 5 is independent.
#
# Afterwards, collapse the scans to one small file to bring home:
#   Rscript scripts_oneoffs/AGE_2024/gather_scans.R

set -euo pipefail

SRC=process/AGE_SY                 # where the SNP list and AGE_SY haplotypes are
DIR=process/AGE_2024               # this experiment
BAMLIST=helpfiles/AGE_Aug13_24/AGE_2024.bams
PARFILE=helpfiles/AGE_Aug13_24/AGE_2024_haplotype_parameters_size75k.R
DESIGN=helpfiles/AGE_Aug13_24/Ageing_Aug13.txt
SMOOTH=100                         # kb -- matches AGE_SY's size75k scans
FOUNDERS=(B1 B2 B3 B4 B5 B6 B7 AB8)

# ── guards ───────────────────────────────────────────────────────────────────
for f in "$BAMLIST" "$PARFILE" "$DESIGN"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done
[ -d "$SRC/Catalog" ] || { echo "ERROR: no $SRC/Catalog to copy the SNP list from" >&2; exit 1; }

missing=0
while read -r b; do
  case "$b" in data/bam/AGE_2024/*) [ -f "$b" ] || { echo "MISSING BAM: $b" >&2; missing=1; } ;; esac
done < "$BAMLIST"
[ "$missing" -eq 0 ] || { echo "aborting -- BAMs not in place" >&2; exit 1; }

# ── 1. SNP list + founder counts ─────────────────────────────────────────────
mkdir -p "$DIR/Calls/counts"
if [ -d "$DIR/Catalog" ]; then
  echo "1. SNP list already at $DIR/Catalog -- leaving it"
else
  cp -a "$SRC/Catalog" "$DIR/Catalog"
  echo "1. copied SNP list: $SRC/Catalog -> $DIR/Catalog"
fi
n_copied=0
for f in "${FOUNDERS[@]}"; do
  s="$SRC/Calls/counts/$f.tsv.gz"; d="$DIR/Calls/counts/$f.tsv.gz"
  if [ -f "$d" ]; then continue; fi
  [ -f "$s" ] || { echo "ERROR: no founder counts at $s" >&2; exit 1; }
  cp -a "$s" "$d"; n_copied=$((n_copied+1))
done
echo "   founder counts in place: $(ls "$DIR"/Calls/counts/*.tsv.gz 2>/dev/null | wc -l) files ($n_copied newly copied)"

# ── 2. call SNPs for the 12 samples ──────────────────────────────────────────
# Count only the samples: the founder columns are the copied counts above.
SAMPLE_BAMS="$DIR/sample_bams.txt"
grep '^data/bam/AGE_2024/' "$BAMLIST" > "$SAMPLE_BAMS"
echo "2. calling SNPs for $(wc -l < "$SAMPLE_BAMS") samples ..."
JID_CALL=$(bash pipeline/scripts/call_samples.sh \
    --catalog "$DIR/Catalog" --bamlist "$SAMPLE_BAMS" --dir "$DIR")
echo "   call_samples merge job: $JID_CALL"

# ── 3. haplotypes ────────────────────────────────────────────────────────────
JID_HAPS=$(bash pipeline/scripts/run_haps.sh --after "$JID_CALL" \
    --parfile "$PARFILE" --dir "$DIR")
echo "3. haplotypes (75 kb / 5 kb) job: $JID_HAPS"

# ── 4. the AGE_2024 scan ─────────────────────────────────────────────────────
out=$(bash pipeline/scripts/run_scan.sh --after "$JID_HAPS" \
    --dir "$DIR" --scan AGE_2024 --design "$DESIGN" --smooth "$SMOOTH")
echo "4. AGE_2024 scan (smooth ${SMOOTH}kb) submitted:"
echo "$out" | sed 's/^/     /'

# ── 5. AGE_SY re-scanned on reps 1-6 ─────────────────────────────────────────
echo "5. AGE_SY reps 1-6 ..."
module load R/4.2.2 2>/dev/null || true
Rscript scripts_oneoffs/AGE_2024/make_reps1_6_designs.R
for scan in AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M; do
  d="helpfiles/AGE_SY/${scan}_R1to6.test.txt"
  o=$(bash pipeline/scripts/run_scan.sh --dir "$SRC" \
        --scan "${scan}_R1to6" --design "$d" --smooth "$SMOOTH")
  echo "   ${scan}_R1to6: $(echo "$o" | tail -1)"
done

cat <<EOF

------------------------------------------------------------------
submitted. steps 2-4 are chained; step 5 runs now.

  process/AGE_2024/Scans/AGE_2024/AGE_2024.scan.txt          <- the pilot
  process/AGE_SY/Scans/AGE_SY{10,20}_{F,M}_R1to6/...         <- matched 6 reps

when everything is done:
  Rscript scripts_oneoffs/AGE_2024/gather_scans.R
  # then scp the one file it names
------------------------------------------------------------------
EOF
