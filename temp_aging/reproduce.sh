#!/bin/bash
#SBATCH --job-name=age_reproduce
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=6G
#SBATCH --time=06:00:00
#SBATCH --output=logs/age_reproduce_%j.out
#
# reproduce.sh — rebuild EVERY derived result from the scans on disk.
#
#   sbatch temp_aging/reproduce.sh        # THIS. it takes minutes, not seconds.
#
# Not `bash` on the head node: the zoom step alone reads four ~257 MB files.
# The SBATCH headers above are already set for it.
#
# RUN ON HPC3 from the repo root. process/ is not on the laptop.
#
# The scans produce two primary things and nothing else is independent:
#
#   Scans/<scan>/<scan>.scan.txt              Wald and h2 per window
#   Scans/<scan>/<scan>.meansBySample.<chr>.txt   founder frequencies per pool
#
# Everything below is a function of those two. So after ANY new scan -- a
# parameter changed, a bug fixed, replicates added -- run this and every number
# and figure input is current. Nothing here resubmits a scan; it derives.
#
# 3 cores for 18 GB: make_zoom_means reads four ~257 MB files. Everything else
# fits in far less. See pipeline/Slurm.md -- standard caps at 6 GB per core.

set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

# Refuse to run on a login node. This reads four ~257 MB files and takes
# minutes; run interactively during the day it earns the lab a warning from
# RCIC. If SLURM_JOB_ID is unset we are not in a job.
if [ -z "${SLURM_JOB_ID:-}" ] && [[ "$(hostname -s)" == login* ]]; then
  cat >&2 <<'MSG'
refusing to run on a login node.

  sbatch temp_aging/reproduce.sh

The SBATCH headers in this file are already set (3 cores, 6 GB each, 6 h).
To override deliberately -- a tiny test, off-hours -- set ALLOW_HEADNODE=1.
MSG
  [ "${ALLOW_HEADNODE:-0}" = "1" ] || exit 1
  echo "ALLOW_HEADNODE=1 set; continuing on $(hostname -s)." >&2
fi

module load R/4.2.2 2>/dev/null || true

SPLITHALF=process/AGE_SY_splithalf
fail=0

step () {   # label, then the command
  local label=$1; shift
  echo
  echo "── $label ──────────────────────────────────────────────"
  if "$@"; then
    echo "   ok"
  else
    echo "   FAILED: $label" >&2
    fail=1
  fi
}

# 1. the two gathered tables: Wald/h2 for the 4 scans, and the 8 split halves
step "gather" Rscript scripts_oneoffs/AGE_SY/nov_only/gather.R

# 2. the variance partition behind Figure 1c. BOTH arguments are required --
#    varcomp_H2.R defaults to the 12-replicate paths and will otherwise rebuild
#    the wrong file, leaving Figure 1c on a stale one with no error.
step "variance partition" \
  Rscript scripts_oneoffs/AGE_SY/splithalf/varcomp_H2.R \
    "${SPLITHALF}/AGE_SY_splithalf_H2_no89.txt.gz" \
    "${SPLITHALF}/H2_varcomp_by_window_no89.txt.gz"

# 3. founder frequencies around the seven zoom peaks, from the means files
step "zoom means" Rscript temp_aging/make_zoom_means.R

# 4. every number quoted in the prose, each with its input provenance
step "numbers" bash temp_aging/run_numbers.sh

echo
echo "════════════════════════════════════════════════════════"
if [ "$fail" -eq 0 ]; then
  echo "all derived results rebuilt."
else
  echo "SOMETHING FAILED -- see above. Do not trust the outputs." >&2
fi
echo
echo "current:"
for f in process/AGE_SY/AGE_SY_4scan_no89.txt.gz \
         "${SPLITHALF}/AGE_SY_splithalf_H2_no89.txt.gz" \
         "${SPLITHALF}/H2_varcomp_by_window_no89.txt.gz" \
         process/AGE_SY/AGE_SY_zoom_means.txt.gz ; do
  if [ -f "$f" ]; then printf '  %-58s %s\n' "$f" "$(date -r "$f" '+%Y-%m-%d %H:%M')"
  else                 printf '  %-58s MISSING\n' "$f"; fi
done
echo
echo "NOT rebuilt here -- the SNP scan is imputed from the smoothed haplotype"
echo "rds, so a new scan makes it stale. Resubmit it separately (no --sex: the"
echo "chrX dosage is already in the Num it reads from the rds):"
echo "  bash scripts_oneoffs/AGE_SY/nov_only/run_snp_scans.sh"
echo
echo "then: git add temp_aging/numbers/ temp_aging/peak_table.txt && git commit && git push"
exit "$fail"
