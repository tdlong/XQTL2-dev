#!/usr/bin/env bash
# pull.sh — bring the small derived files down from HPC3.
#
#   bash temp_aging/pull.sh
#
# Run on the laptop after a scan or a reproduce run. Same command every time.
# One rsync, ~29 MB, everything the figures and the numbers scripts read except
# the SNP scan, which is 81 MB and feeds Figure 3 alone -- add --snp for that.
#
# With these local, figures redraw in seconds:  bash temp_aging/make_figures.sh --no-fetch

set -uo pipefail
[ -d scripts_oneoffs/AGE_SY ] || { echo "run from the XQTL2-dev root" >&2; exit 1; }

REMOTE=tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/process/
mkdir -p process/AGE_SY/Calls process/AGE_SY_splithalf

{
  cat <<'SMALL'
AGE_SY/AGE_SY_4scan_no89.txt.gz
AGE_SY/AGE_SY_zoom_means.txt.gz
AGE_SY/Calls/refalt_qc.txt
AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz
AGE_SY_splithalf/H2_varcomp_by_window_no89.txt.gz
SMALL
  [ "${1:-}" = "--snp" ] && echo AGE_SY/AGE_SY_4snpscan_no89.txt.gz
} | rsync -av --files-from=- "$REMOTE" process/

echo
du -sh process/* 2>/dev/null
