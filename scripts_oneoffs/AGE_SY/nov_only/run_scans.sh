#!/bin/bash
# run_scans.sh — the twelve AGE_SY scans without replicates 8 and 9, and
# everything downstream of them.
#
# ONE command, from the repo root on HPC3, after `git pull`:
#   bash scripts_oneoffs/AGE_SY/nov_only/run_scans.sh
#
# It submits the twelve scans, then chains gather, the zoom means and
# run_numbers on them with SLURM dependencies and returns. There is nothing to
# run by hand afterwards and nothing to wait around for -- when the last job
# clears, temp_aging/numbers/ and the three process/ files are current.
#
# Replicates 8 and 9 are the May 2023 cage; dropping them leaves a single-cage
# experiment on 10 replicates, splitting evenly 5/5 into odd and even.
#
#   4 x <scan>_no89              -> process/AGE_SY/Scans            (10 reps)
#   8 x <scan>_no89_{odd,even}   -> process/AGE_SY_splithalf/Scans  (5 reps each)
#
# Haplotypes already exist in both folders, so this is scan-only. The split-half
# folder reads its Haps through a symlink to process/AGE_SY, set up by
# scripts_oneoffs/AGE_SY/splithalf/run_splithalf_scans.sh -- it must already be
# there, and this script checks.
#

set -euo pipefail

SMOOTH=100                       # kb -- matches every other AGE_SY scan

# Every AGE_SY pool is one cage x treatment x sex, so each scan is single sex and
# takes --sex from the last field of its name. Without it run_scan.sh defaults to
# "mixed" (0.75 X chromosomes per fly), which is right for a 50:50 pool and wrong
# for these: 1.0 for females, 0.5 for males (XQTL2 #38). chrX only.
FULL_DIR=process/AGE_SY
HALF_DIR=process/AGE_SY_splithalf
DESIGNS=helpfiles/AGE_SY/nov_only
SCANS=(AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M)

jobid_from() {
  local id
  id=$(printf '%s\n' "$1" | grep -oE '[0-9]+' | tail -1)
  [[ -n "$id" ]] || { echo "ERROR: no job id in:" >&2
                      printf '%s\n' "$1" >&2; return 1; }
  printf '%s' "$id"
}

[ -d "$FULL_DIR/Haps" ] || { echo "ERROR: no haplotypes at $FULL_DIR/Haps" >&2; exit 1; }
[ -e "$HALF_DIR/Haps" ] || { echo "ERROR: no Haps symlink at $HALF_DIR/Haps" >&2; exit 1; }

echo "building designs ..."
module load R/4.2.2 2>/dev/null || true
Rscript scripts_oneoffs/AGE_SY/nov_only/make_designs.R

echo
echo "submitting ..."
CONCAT_IDS=()

submit() {   # name, dir, design, sex
  [ -f "$3" ] || { echo "ERROR: missing design $3" >&2; exit 1; }
  local out id; out=$(bash pipeline/scripts/run_scan.sh \
      --dir "$2" --scan "$1" --design "$3" --smooth "$SMOOTH" --sex "$4")
  id=$(jobid_from "$out")
  CONCAT_IDS+=("$id")
  echo "   $(printf '%-26s' "$1") sex $4  concat $id"
}

for s in "${SCANS[@]}"; do
  submit "${s}_no89" "$FULL_DIR" "$DESIGNS/${s}.no89.txt" "${s##*_}"
done
for s in "${SCANS[@]}"; do
  for h in odd even; do
    submit "${s}_no89_${h}" "$HALF_DIR" "$DESIGNS/${s}.no89.${h}.txt" "${s##*_}"
  done
done

# ---------------------------------------------------------------------------
# Everything downstream, chained on the twelve concat jobs. Nothing below needs
# a human between steps, so it is submitted now and runs when the scans land.
# ---------------------------------------------------------------------------
DEP_SCANS="afterok:$(IFS=:; echo "${CONCAT_IDS[*]}")"
SB="sbatch --parsable -A tdlong_lab -p standard"

# Everything derived is one script, so the definition lives in one place and
# `bash temp_aging/reproduce.sh` rebuilds the same things without resubmitting
# scans. 3 cores for 18 GB -- the zoom step reads four ~257 MB files.
JID_DERIVE=$($SB --dependency="$DEP_SCANS" \
    --job-name=age_reproduce --cpus-per-task=3 --mem-per-cpu=6G --time=06:00:00 \
    --output=logs/age_reproduce_%j.out \
    --wrap="bash temp_aging/reproduce.sh")

cat <<EOF

------------------------------------------------------------------
submitted, chained end to end. nothing else to run.

  12 scans      ${CONCAT_IDS[0]} .. ${CONCAT_IDS[${#CONCAT_IDS[@]}-1]}
  derive        $JID_DERIVE   (after the scans: gather, partition, zoom, numbers)

watch:    squeue -u \$USER
log:      logs/age_reproduce_<jobid>.out

when $JID_DERIVE finishes, temp_aging/numbers/ is current and these exist:
  process/AGE_SY/AGE_SY_4scan_no89.txt.gz
  process/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz
  process/AGE_SY/AGE_SY_zoom_means.txt.gz
  process/AGE_SY_splithalf/H2_varcomp_by_window_no89.txt.gz

STILL TO RUN BY HAND -- the SNP scans are NOT chained here. They read the
smoothed haplotypes these scans just rewrote, so they are now stale:

  bash scripts_oneoffs/AGE_SY/nov_only/run_snp_scans.sh

then, to bring it back:
  git add temp_aging/numbers/ && git commit -m "numbers after the chrX sex fix" && git push
------------------------------------------------------------------
EOF
