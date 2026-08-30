#!/usr/bin/env bash
# make_figures.sh — fetch what the figures need and draw all of them.
#
#   bash temp_aging/make_figures.sh            # fetch, then draw everything
#   bash temp_aging/make_figures.sh --no-fetch # draw from what is already here
#   bash temp_aging/make_figures.sh 1 3        # only figures 1 and 3
#
# RUN ON THE LAPTOP, from the repo root, on VPN. One rsync, one connection --
# not a pile of scp lines. Existing local copies are overwritten: after a rerun
# on the cluster they are stale, which is the whole point of fetching.

set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

REMOTE=tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/process/
FETCH=1; WANT=()
for a in "$@"; do
  case $a in
    --no-fetch) FETCH=0 ;;
    1|2|3|rr)   WANT+=("$a") ;;
    *) echo "usage: make_figures.sh [--no-fetch] [1] [2] [3] [rr]" >&2; exit 1 ;;
  esac
done
[ ${#WANT[@]} -eq 0 ] && WANT=(1 2 3 rr)

if [ "$FETCH" -eq 1 ]; then
  echo "fetching from HPC3 ..."
  mkdir -p process/AGE_SY/Calls process/AGE_SY_splithalf
  rsync -av --files-from=- "$REMOTE" process/ <<'LIST'
AGE_SY/AGE_SY_4scan_no89.txt.gz
AGE_SY/AGE_SY_zoom_means.txt.gz
AGE_SY/AGE_SY_4snpscan_no89.txt.gz
AGE_SY/Calls/refalt_qc.txt
AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz
AGE_SY_splithalf/H2_varcomp_by_window_no89.txt.gz
LIST
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "rsync failed (rc=$rc). On VPN? Missing files are listed above." >&2
    echo "Carrying on with whatever is already local." >&2
  fi
fi

echo
for n in "${WANT[@]}"; do
  case $n in
    rr) script=temp_aging/make_figure_rr.R ;;
    *)  script=temp_aging/make_figure${n}.R ;;
  esac
  echo "── figure $n ─────────────────────────────────"
  # Figure 1 is the genetic axis. make_figure1.R defaults to X_UNIT=Mb, but the
  # paper uses cM -- Figure1_cM_plot.png, linkage groups concatenated. Set it
  # here so the default run produces the figure that is actually used.
  if [ "$n" = "1" ]; then export X_UNIT=cM; else unset X_UNIT; fi
  if Rscript "$script"; then echo "   ok"; else echo "   FAILED: $script" >&2; fail=1; fi
done

echo
echo "written:"
ls -lt temp_aging/*.png 2>/dev/null | awk '{printf "  %-34s %s %s %s\n", $9, $6, $7, $8}' | head
exit "${fail:-0}"
