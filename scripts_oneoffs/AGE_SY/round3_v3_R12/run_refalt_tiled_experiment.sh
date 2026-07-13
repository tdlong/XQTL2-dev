#!/bin/bash
###############################################################################
# run_refalt_tiled_experiment.sh — AGE_SY round 3 (v3): validate the new caller
#
# EXPERIMENT: prove the tiled/scatter SNP caller (bam2bcf2REFALT.tiled.sh)
# produces byte-identical RefAlt.<chr>.txt to the validated caller.
#
#   Reference (validated): process/AGE_SY_v3       <- job 54049866 (run_refalt.sh)
#   Candidate  (new/tiled): process/AGE_SY_v3_tiled <- submitted here
#
# Both use the SAME bam list (helpfiles/AGE_SY/AGE_SY.bams) so inputs are
# identical. The candidate writes to a SEPARATE dir and only reads the bams,
# so it runs safely alongside the still-running reference job.
#
# This submits ~29 tile tasks + a dependent reassembly job, then prints the
# reassembly job ID. When BOTH that job and the reference (54049866) have
# finished, run the compare command printed at the end.
#
# Run from: XQTL2-dev repo root, on the cluster.
###############################################################################
set -euo pipefail

BAMLIST=helpfiles/AGE_SY/AGE_SY.bams
REF_DIR=process/AGE_SY_v3
CAND_DIR=process/AGE_SY_v3_tiled

[ -f "$BAMLIST" ] || { echo "ERROR: $BAMLIST not found (run from repo root, on the cluster)" >&2; exit 1; }

# --- submit the candidate (tiled) caller into its own directory ---
JID=$(bash pipeline/scripts/run_refalt.tiled.sh \
        --bamlist "$BAMLIST" \
        --dir     "$CAND_DIR")
echo "candidate (tiled) reassembly job: $JID   -> $CAND_DIR"

cat <<EOF

================================================================
Candidate submitted. Monitor:  squeue -u \$USER
  (waits on ~29 tile tasks, then reassembles into $CAND_DIR/RefAlt.<chr>.txt)

Reference run: job 54049866 -> $REF_DIR   (let it finish)

--- Compare ONCE BOTH job $JID AND job 54049866 are done: ---
  bash pipeline/scripts/compare_refalt.sh $REF_DIR $CAND_DIR

PASS = IDENTICAL on every chromosome.
Anything DIFFERS prints the first offending lines -> a real bug to chase.
================================================================
EOF
