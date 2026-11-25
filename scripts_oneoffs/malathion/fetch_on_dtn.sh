#!/bin/bash
###############################################################################
# fetch_on_dtn.sh — start the malathion download on HPC3's data-transfer host
# from the login node, without logging in there yourself.
#
# Transfers run on access-hpc3.rcic.uci.edu and nowhere else: not on login nodes,
# not on compute nodes, not in Slurm (https://rcic.uci.edu/storage/data-transfer.html,
# and pipeline/Slurm.md, "Downloads"). That host mounts the same filesystem, so
# this ssh's there, cd's into this same repo path, and starts fetch_sra.sh
# detached with nohup. ssh may ask you to authenticate once. Then it returns.
#
# Watch from the login node (the log is on the shared filesystem):
#   tail logs/malathion/fetch_dtn.out
# It ends with "done: 16 / 16 files verified, exit 0". If fewer, run this again:
# only missing or corrupt files are refetched. Never start it twice at once.
#
# Run from the repo root on a login node:
#   bash scripts_oneoffs/malathion/fetch_on_dtn.sh
###############################################################################
set -euo pipefail
[ -d scripts_oneoffs/malathion ] || { echo "ERROR: run from the XQTL2-dev root" >&2; exit 1; }
DTN=access-hpc3.rcic.uci.edu
REPO=$PWD          # same path on the DTN: DFS is mounted there too
LOG=logs/malathion/fetch_dtn.out
mkdir -p logs/malathion

if pgrep -u "$USER" -f fetch_sra.sh >/dev/null 2>&1; then
  echo "ERROR: fetch_sra.sh is already running under your user here; not starting another" >&2; exit 1
fi
if [ -f "$LOG" ] && tail -1 "$LOG" | grep -q "verified, exit 0"; then
  n_ok=$(ls data/raw/malathion_2022/*.md5ok 2>/dev/null | wc -l)
  [ "$n_ok" -eq 16 ] && { echo "already done: 16 / 16 verified (see $LOG)"; exit 0; }
fi

echo "starting fetch_sra.sh on $DTN in $REPO (you may be asked to authenticate)"
ssh "$DTN" "cd '$REPO' && nohup bash scripts_oneoffs/malathion/fetch_sra.sh </dev/null >/dev/null 2>&1 &"
echo "started. watch it with:   tail $LOG"
echo "when it reports 16 / 16 verified:   bash scripts_oneoffs/malathion/rerun_2022_from_reads.sh"
