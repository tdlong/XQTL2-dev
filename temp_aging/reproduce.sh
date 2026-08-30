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
# Run from the XQTL2-dev root, like everything else in this repo. Scripts live
# in XQTL2 and are reached through the pipeline symlink; paths below are all
# relative to the root. No cd -- sbatch starts you where you submitted from.
if [ ! -d scripts_oneoffs/AGE_SY ]; then
  echo "run this from the XQTL2-dev root; cwd is $(pwd)" >&2; exit 1
fi

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

# SCOPE main rebuilds only what the four main scans support: the gathered Wald
# and h2 table, the zoom means, and the corrected h2. It deliberately skips the
# variance partition and run_numbers, because both read the split-half file --
# and pairing freshly rerun main scans with split halves from an older run is
# exactly the kind of silent mixing that has bitten this project.
# h2: only re-derive the heritability from the existing smoothed rds. Nothing
# upstream changes when the h2 estimator does, so a rescan is not needed.
SCOPE=${1:-main}
case $SCOPE in main|all|h2) ;; *) echo "usage: reproduce.sh [main|all|h2]" >&2; exit 1 ;; esac
echo "scope:    $SCOPE"
echo "repo:     $(git log -1 --format=%h) $(git log -1 --format=%s | cut -c1-52)"
echo "pipeline: $(git -C pipeline log -1 --format=%h 2>/dev/null) $(git -C pipeline log -1 --format=%s 2>/dev/null | cut -c1-52)"

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
[ "$SCOPE" = h2 ] || step "gather" Rscript scripts_oneoffs/AGE_SY/nov_only/gather.R

# 2. the variance partition behind Figure 1c. BOTH arguments are required --
#    varcomp_H2.R defaults to the 12-replicate paths and will otherwise rebuild
#    the wrong file, leaving Figure 1c on a stale one with no error.
if [ "$SCOPE" = all ]; then
step "variance partition" \
  Rscript scripts_oneoffs/AGE_SY/splithalf/varcomp_H2.R \
    "${SPLITHALF}/AGE_SY_splithalf_H2_no89.txt.gz" \
    "${SPLITHALF}/H2_varcomp_by_window_no89.txt.gz"
else
  echo; echo "── variance partition ── skipped (SCOPE=main; needs the split halves)"
fi

# 3. founder frequencies around the seven zoom peaks, from the means files
[ "$SCOPE" = h2 ] || step "zoom means" Rscript temp_aging/make_zoom_means.R

# 4. bias-corrected Falconer h2, one file per scan. The pipeline's own h2 columns
#    are the uncorrected Cutler ones, whose bias term saturates against the
#    penetrance clamp and so under-corrects exactly where the variance is largest
#    -- the male X (XQTL2 #40). h2_from_scan.R subtracts the exact Falconer bias
#    instead. Needs the smoothed rds and the design, so it runs here, not on
#    meansBySample.
for sc in AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M; do
  step "h2: ${sc}_no89" \
    Rscript pipeline/scripts/h2_from_scan.R \
      --dir   process/AGE_SY \
      --scan  "${sc}_no89" \
      --rfile "helpfiles/AGE_SY/nov_only/${sc}.no89.txt"
done

# The eight half-scans, which the variance partition consumes. h2_rep must be
# computed per half on 5 replicates -- the partition cannot take the 10-replicate
# value and split it. Skipped under SCOPE=main, where the halves are not rerun.
if [ "$SCOPE" != main ]; then
  for sc in AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M; do
    for hf in odd even; do
      step "h2: ${sc}_no89_${hf}" \
        Rscript pipeline/scripts/h2_from_scan.R \
          --dir   process/AGE_SY_splithalf \
          --scan  "${sc}_no89_${hf}" \
          --rfile "helpfiles/AGE_SY/nov_only/${sc}.no89.${hf}.txt"
    done
  done
fi

# 5. the R2 smoothing correction must exist before the h2 above is trusted:
#    without <scan>.smooth_r2.txt, hap_scan silently applies R2=1 and the h2
#    bias over-subtracts by 1/R2. run_scan.sh submits the diagnostic, so a
#    missing file means that job failed rather than that it was skipped.
echo
echo "── R2 smoothing correction ──────────────────────────────────"
missing_r2=0
for sc in AGE_SY10_F AGE_SY20_F AGE_SY10_M AGE_SY20_M; do
  f="process/AGE_SY/Scans/${sc}_no89/${sc}_no89.smooth_r2.txt"
  if [ -f "$f" ]; then printf '   %-18s R2 = %s\n' "$sc" "$(cat "$f")"
  else                 printf '   %-18s MISSING\n' "$sc"; missing_r2=1; fi
done
if [ "$missing_r2" -eq 1 ]; then
  echo "   ^ hap_scan applied no correction for these; the h2 bias is too large" >&2
  echo "     by 1/R2. Check the smooth_r2 job in the run_scan.sh chain." >&2
  fail=1
fi

# 5b. can the partition run on h2_rep? Prints; run_numbers.sh captures it to
#     numbers/partition_check.txt, which comes back through git. The point of
#     doing it here is that the per-window per-half h2 is far too large to move,
#     while the answer is a few dozen lines.
# Written straight to numbers/, not just printed: under scope h2 run_numbers.sh
# does not run, so relying on it to capture this meant the answer stayed in the
# job log and never came back through git.
if [ "$SCOPE" != main ]; then
  mkdir -p temp_aging/numbers
  {
    echo "# partition_check.R"
    echo "# run:      $(date '+%Y-%m-%d %H:%M')"
    echo "# script:   $(git log -1 --format=%h -- temp_aging/partition_check.R 2>/dev/null || echo uncommitted)"
    echo "# pipeline: $(git -C pipeline log -1 --format=%h 2>/dev/null)"
    echo "#"
  } > temp_aging/numbers/partition_check.txt
  step "partition check" bash -c \
    'Rscript temp_aging/partition_check.R >> temp_aging/numbers/partition_check.txt 2>&1'
  sed -n '7,80p' temp_aging/numbers/partition_check.txt
fi

# 6. every number quoted in the prose, each with its input provenance. Skipped
#    under SCOPE=main: significant_regions.R reads the split-half file, so the
#    numbers would mix new main scans with stale halves.
if [ "$SCOPE" = all ]; then
  step "numbers" bash temp_aging/run_numbers.sh
else
  echo; echo "── numbers ── skipped (SCOPE=main; several read the split halves)"
fi

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
