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
# Where the repo is. Under sbatch, $0 is SLURM's spool copy of this script, not
# this file, so dirname "$0"/.. points somewhere unrelated and every Rscript
# below fails with "cannot open file". SLURM_SUBMIT_DIR is the directory sbatch
# was run from, which is the repo root.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  cd "$SLURM_SUBMIT_DIR" || exit 1
else
  cd "$(dirname "$0")/.." || exit 1
fi

# Prove it, however we got here.
if [ ! -f temp_aging/reproduce.sh ] || [ ! -d scripts_oneoffs/AGE_SY ]; then
  echo "not at the repo root -- cwd is $(pwd)" >&2
  echo "sbatch this from the repo root: cd /dfs7/adl/tdlong/fly_pool/XQTL2-dev" >&2
  exit 1
fi
echo "repo: $(pwd)"

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

# Load R and prove it worked. In a batch job `module` may not be initialised the
# way it is in a login shell, and hiding that behind `2>/dev/null || true` makes
# every step below fail identically with no clue why.
if ! command -v module >/dev/null 2>&1; then
  for init in /opt/rcic/Modules/default/init/bash /usr/share/Modules/init/bash \
              /etc/profile.d/modules.sh; do
    [ -r "$init" ] && . "$init" && break
  done
fi
module load R/4.2.2 || echo "module load R/4.2.2 failed" >&2

if ! command -v Rscript >/dev/null 2>&1; then
  cat >&2 <<MSG
Rscript is not on PATH, so every step below would fail for that reason alone.

  module list          # what is loaded
  module avail R       # what R versions exist
  which Rscript

PATH=$PATH
MSG
  exit 1
fi
echo "using $(command -v Rscript)  [$(Rscript -e 'cat(R.version.string)' 2>&1 | head -1)]"

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

# 4. bias-corrected Falconer h2, one file per scan. The pipeline's own h2 columns
#    are the uncorrected Cutler ones, whose bias term saturates against the
#    penetrance clamp and so under-corrects exactly where the variance is largest
#    -- the male X (XQTL2 #40). h2_from_scan.R subtracts the exact Falconer bias
#    instead. Needs the smoothed rds and the design, so it runs here, not on
#    meansBySample.
for sc in AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M; do
  step "falconer h2: ${sc}_no89" \
    Rscript pipeline/scripts/h2_from_scan.R \
      --dir   process/AGE_SY \
      --scan  "${sc}_no89" \
      --rfile "helpfiles/AGE_SY/nov_only/${sc}.no89.txt"
done

# 5. every number quoted in the prose, each with its input provenance
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
         process/AGE_SY/AGE_SY_zoom_means.txt.gz \
         process/AGE_SY/Scans/AGE_SY20_M_no89/AGE_SY20_M_no89.h2_falconer.txt ; do
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
