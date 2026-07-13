#!/bin/bash
#SBATCH --job-name=compare_refalt
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2:00:00
#SBATCH -o compare_refalt_experiment.%j.out
###############################################################################
# compare_refalt_experiment.sh — AGE_SY r3: does the tiled caller match?
#
# Runs pipeline/scripts/compare_refalt.sh on the two experiment output dirs and
# writes the verdict to compare_refalt_experiment.<jobid>.out — one line per
# chromosome: IDENTICAL / REORDERED / DIFFERS / MISSING.
#
#   reference (validated): process/AGE_SY_v3        (job 54049866)
#   candidate (tiled):     process/AGE_SY_v3_tiled   (tiled reassembly job)
#
# WHEN TO SUBMIT: only after BOTH runs have finished. Check `squeue -u $USER`
# — both the reference and the tiled reassembly job must be gone. Submitting
# early just reports MISSING for the tiled files (harmless, but pointless).
#
# Runs as a batch job (not on the login node) because compare_refalt.sh sorts
# the full RefAlt tables for any chromosome that is not byte-identical.
#
# Submit from XQTL2-dev repo root, on the cluster:
#   sbatch scripts_oneoffs/AGE_SY/round3_v3_R12/compare_refalt_experiment.sh
###############################################################################
set -uo pipefail

REF_DIR=process/AGE_SY_v3
CAND_DIR=process/AGE_SY_v3_tiled

echo "Comparing RefAlt tables:"
echo "  reference (validated): $REF_DIR"
echo "  candidate (tiled):     $CAND_DIR"
echo

bash pipeline/scripts/compare_refalt.sh "$REF_DIR" "$CAND_DIR"
rc=$?

echo
if [ "$rc" -eq 0 ]; then
    echo "RESULT: PASS — every chromosome IDENTICAL (REORDERED counts as a warning, not a fail)."
else
    echo "RESULT: FAIL — a chromosome DIFFERS or is MISSING (see the lines above)."
fi
exit $rc
