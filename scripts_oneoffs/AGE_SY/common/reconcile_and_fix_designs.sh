#!/usr/bin/env bash
# reconcile_and_fix_designs.sh
#
# Purpose: unblock the pull and CORRECT THE DESIGN FILES.
#
# The cluster tree is dirty, so `git pull` refuses. The dirty files split into:
#   - AGE_SY.bams + AGE_SY_haplotype_parameters.R  = REAL R12 work -> COMMIT
#   - AGE_SY*.test.txt (4 design files)            = DERIVED from the xlsx,
#                                                    about to be regenerated -> DISCARD
# Then pull (gets the new hard-fail generator), regenerate the 4 design files
# from summary_info_v1.xlsx (which has R12), and commit them locally.
#
# Does NOT push. Does NOT touch AGE_SY.bams contents on disk (a running SNP job
# read it at launch; git commit records, it does not rewrite the file bytes).
# Aborts on the first surprise instead of plowing ahead.
#
# Run on the cluster from repo root:
#   bash scripts_oneoffs/AGE_SY/common/reconcile_and_fix_designs.sh 2>&1 | tee ~/fix.txt

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
say(){ printf '\n>>> %s\n' "$*"; }

REAL=(helpfiles/AGE_SY/AGE_SY.bams helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R)
DESIGN=(helpfiles/AGE_SY/AGE_SY10_F.test.txt helpfiles/AGE_SY/AGE_SY10_M.test.txt \
        helpfiles/AGE_SY/AGE_SY20_F.test.txt helpfiles/AGE_SY/AGE_SY20_M.test.txt)
XLSX=helpfiles/AGE_SY/summary_info_v1.xlsx
GEN=scripts_oneoffs/AGE_SY/common/make_AGE_SY_design_files.py

# --- precondition: branch + clean-ish expectations -------------------------
[ "$(git rev-parse --abbrev-ref HEAD)" = "main" ] || { echo "ABORT: not on main"; exit 1; }

say "Precondition: xlsx (data source) must be UNMODIFIED"
git diff --quiet -- "$XLSX" || { echo "ABORT: $XLSX is modified — inspect before proceeding"; exit 1; }
echo "  xlsx clean (cluster == committed). OK"

say "Precondition: the two 'real work' files must actually be modified"
for f in "${REAL[@]}"; do
  git diff --quiet -- "$f" && { echo "ABORT: expected $f to be modified but it is clean"; exit 1; }
done
echo "  AGE_SY.bams + haplotype_parameters.R are modified as expected. OK"

# --- step 1: discard the derived design files ------------------------------
say "Step 1: discard the 4 derived design files (they get regenerated below)"
git checkout -- "${DESIGN[@]}" 2>/dev/null || true   # if some are already clean, fine
echo "  reverted: ${DESIGN[*]}"

# --- step 2: commit the real R12 work --------------------------------------
say "Step 2: commit AGE_SY.bams + names_in_bam (R12 work from 03_update_helpfiles.sh)"
git add "${REAL[@]}"
git commit -q -m "AGE_SY R12: add R12 BAMs to AGE_SY.bams + names_in_bam"
echo "  committed $(git rev-parse --short HEAD)"

# --- step 3: pull the new (hard-fail) generator ----------------------------
say "Step 3: rebase onto origin/main (brings the hard-fail generator)"
git fetch origin main
git rebase origin/main
echo "  now at $(git rev-parse --short HEAD); behind=$(git rev-list --count HEAD..origin/main)"

# --- step 4: regenerate the design files from the xlsx ---------------------
say "Step 4: regenerate the 4 design files (hard-fails if any R12 count is missing)"
python3 "$GEN"

# --- step 5: commit the corrected design files -----------------------------
say "Step 5: commit the corrected design files"
git add "${DESIGN[@]}"
if git diff --cached --quiet; then
  echo "  (no change vs committed design files — nothing to commit)"
else
  git commit -q -m "AGE_SY R12: regenerate design files from summary_info_v1.xlsx"
  echo "  committed $(git rev-parse --short HEAD)"
fi

# --- report ----------------------------------------------------------------
say "RESULT"
git --no-pager log --oneline -4
echo "--- status ---"
git status --short
echo "--- R12 rows present in each regenerated design file ---"
for f in "${DESIGN[@]}"; do
  printf '  %s : %s\n' "$f" "$(grep -c '"R12"' "$f" || echo 0) R12 row(s)"
done
echo
echo "NOT pushed yet. Review above, then push with:"
echo "  git push origin main"
