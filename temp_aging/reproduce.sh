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
SCOPE=${1:-main}
case $SCOPE in main|all) ;; *) echo "usage: reproduce.sh [main|all]" >&2; exit 1 ;; esac
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
step "gather" Rscript scripts_oneoffs/AGE_SY/nov_only/gather.R

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
step "zoom means" Rscript temp_aging/make_zoom_means.R

# H2 and H2_vc now come straight out of hap_scan.R in scan.txt (XQTL2 b93b98b);
# h2_from_scan.R was deleted upstream, so there is no separate h2 step here any
# more. gather.R picks the columns up.

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

# 5c. the figures. Drawn here rather than on the laptop: they read process/,
#     which is only here, and they are committed so they come back through git
#     with the numbers instead of over a connection that keeps failing.
if [ "$SCOPE" != main ]; then
  for fg in 1 2 3 rr; do
    case $fg in rr) scr=temp_aging/make_figure_rr.R ;; *) scr=temp_aging/make_figure${fg}.R ;; esac
    if [ "$fg" = 1 ]; then export X_UNIT=cM; else unset X_UNIT; fi
    step "figure ${fg}" Rscript "$scr"
  done
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
         ; do
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
