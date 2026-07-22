#!/usr/bin/env bash
# cluster_sync.sh — after ANY cluster activity, run this ONE command (always
# identical, no arguments):
#
#     bash scripts/cluster_sync.sh
#
# Claude Code runs off-cluster and sees only git. This captures everything it
# needs about a run and pushes it, so Claude `git pull`s and reads it — no pasting.
# It captures, generically (it does not need to know which project you ran):
#   - all pipeline job logs (slurm-<jobid>.out dropped in the repo root)
#   - SLURM job state (sacct history + squeue)
#   - a manifest (names/sizes) of everything under process/
# All into the tracked logs/ tree, then commit + push.
#
# Run from the repo root on the cluster. Everything under logs/ is tracked text;
# delete old files anytime.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
mkdir -p logs
STAMP=$(date '+%Y%m%d_%H%M%S')

# 1. sweep pipeline default logs (slurm-<jobid>.out) from the repo root into logs/
shopt -s nullglob
n=0; for f in slurm-*.out; do mv "$f" "logs/"; n=$((n+1)); done
[[ $n -gt 0 ]] && echo "swept $n slurm-*.out into logs/"

# 2. SLURM job state: recent history + current queue
{
  echo "### sacct (since yesterday)  $STAMP"
  sacct -X --starttime "$(date -d 'yesterday' +%F 2>/dev/null || date -v-1d +%F)" \
     --format=JobID,JobName%20,State,ExitCode,Elapsed,Start,End 2>/dev/null
  echo; echo "### squeue -u $USER"
  squeue -u "$USER" 2>/dev/null
} > "logs/state_$STAMP.txt"

# 3. manifest of process/ output (names + sizes, not the bulk data)
{
  for d in process/*/; do
    [[ -d "$d" ]] || continue
    echo "### $d"; ls -la "$d" 2>/dev/null | head -80; echo
  done
} > "logs/manifest_$STAMP.txt"
echo "captured state + manifest ($STAMP)"

# 4. commit + push (rebase first: the Mac side also pushes to this origin)
git add logs/
if git diff --cached --quiet; then echo "nothing new to sync"; exit 0; fi
git commit -q -m "cluster sync $STAMP"
git pull --rebase -q origin main 2>/dev/null || true
git push -q origin main && echo "synced -> git.  Claude: git pull && read logs/"
