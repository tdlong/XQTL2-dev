#!/bin/bash
# run_snp_scans.sh — the four SNP scans for the 10-replicate dataset.
#
# Run from the repo root on HPC3, after `git pull` in BOTH repos:
#   cd ../XQTL2 && git pull && cd -      # the SNP scan changed, see XQTL2 #35
#   bash scripts_oneoffs/AGE_SY/nov_only/run_snp_scans.sh
#
# Imputes per-SNP ALT frequency from the smoothed haplotype estimates and runs a
# df=1 Wald test at every SNP:
#
#     f_ALT(pool, SNP) = h(pool, window) . s(SNP)
#
# Both halves come from Calls/RefAlt.<chr>.txt. Since XQTL2 #35 there is no SNP
# table to prepare and no --founders flag: founder states are read straight from
# the founder columns of RefAlt, and founder names come from --parfile, the same
# file REFALT2haps used. That is also why the founders must have been counted
# into RefAlt in Step 3 -- this script checks before submitting anything.
#
# Uses the same --scan names as the haplotype scans, so output lands beside them:
#   process/AGE_SY/Scans/<scan>_no89/<scan>_no89.snp_scan.txt
#
# The gather is chained on all four scans, so there is nothing to run afterwards.
# It writes process/AGE_SY/AGE_SY_4snpscan_no89.txt.gz and the script prints the
# scp command to fetch it.
#
# The scan is NOT independent of the haplotype scan -- same smoothed haplotypes.

set -uo pipefail

DIR=process/AGE_SY
PARFILE=helpfiles/AGE_SY/AGE_SY_haplotype_parameters_size75k.R
DESIGNS=helpfiles/AGE_SY/nov_only
SCANS=(AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M)
FOUNDERS=(B1 B2 B3 B4 B5 B6 B7 AB8)

jobid_from() {
  local id
  id=$(printf '%s\n' "$1" | grep -oE '[0-9]+' | tail -1)
  [[ -n "$id" ]] || { echo "ERROR: no job id in:" >&2
                      printf '%s\n' "$1" >&2; return 1; }
  printf '%s' "$id"
}

# ── guards ───────────────────────────────────────────────────────────────────
[ -f "$PARFILE" ] || { echo "ERROR: missing $PARFILE" >&2; exit 1; }
for s in "${SCANS[@]}"; do
  d="$DESIGNS/${s}.no89.txt"
  [ -f "$d" ] || { echo "ERROR: missing $d -- run make_designs.R" >&2; exit 1; }
  h="$DIR/Scans/${s}_no89"
  [ -d "$h" ] || { echo "ERROR: no haplotype scan at $h -- run run_scans.sh first" >&2; exit 1; }
done

# The founders have to be IN RefAlt, or there is nothing to read s(SNP) from.
REFALT="$DIR/Calls/RefAlt.chr2L.txt"
[ -f "$REFALT" ] || { echo "ERROR: missing $REFALT" >&2; exit 1; }
hdr=$(head -1 "$REFALT")
missing=()
for f in "${FOUNDERS[@]}"; do
  case "$hdr" in *"REF_$f"*) ;; *) missing+=("$f") ;; esac
done
if [ ${#missing[@]} -gt 0 ]; then
  cat >&2 <<EOF
ERROR: RefAlt has no founder columns for: ${missing[*]}

The SNP scan reads founder states from RefAlt, so the founder BAMs must have been
part of the --bamlist given to call_samples.sh. To add them without recounting
the samples, rerun call_samples.sh with a bamlist of JUST the founders -- their
counts land beside the existing ones and RefAlt is re-merged:

  bash pipeline/scripts/call_samples.sh \\
      --catalog $DIR/Catalog \\
      --bamlist <file listing only the 8 founder BAMs> \\
      --dir     $DIR
EOF
  exit 1
fi
echo "RefAlt carries all 8 founder columns."

# ── submit ───────────────────────────────────────────────────────────────────
echo "submitting ..."
JIDS=()
for s in "${SCANS[@]}"; do
  out=$(bash pipeline/scripts/run_snp_scan.sh \
          --design  "$DESIGNS/${s}.no89.txt" \
          --dir     "$DIR" \
          --scan    "${s}_no89" \
          --parfile "$PARFILE")
  jid=$(jobid_from "$out")
  JIDS+=("$jid")
  echo "   $(printf '%-22s' "${s}_no89") $jid"
done

# ── the gather, chained on all four ──────────────────────────────────────────
mkdir -p logs/AGE_SY
DEP=$(IFS=:; echo "${JIDS[*]}")
GID=$(sbatch --parsable --dependency=afterok:"$DEP" \
        -A tdlong_lab -p standard --time=1:00:00 --mem-per-cpu=12G \
        -J snpgather -o logs/AGE_SY/snp_gather.out \
        --wrap="module load R/4.2.2; Rscript scripts_oneoffs/AGE_SY/nov_only/gather_snp.R")
echo "   $(printf '%-22s' gather) $GID   (afterok on the four above)"

cat <<EOF

------------------------------------------------------------------
4 SNP scans + the gather submitted. Nothing to run afterwards.

When job $GID has finished, from your LAPTOP in the XQTL2-dev repo root:

  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/$DIR/AGE_SY_4snpscan_no89.txt.gz $DIR/
  Rscript temp_aging/make_figure3.R

Gather log (row count and file size): logs/AGE_SY/snp_gather.out
------------------------------------------------------------------
EOF
