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

# 3. manifest of process/ output: a REAL file inventory (type, size, mtime, path) 4
#    levels deep, so Claude KNOWS what completed where — which catalogs exist,
#    whether Calls/RefAlt and Calls/counts/* are there — without anyone ls-ing by hand.
find process -maxdepth 4 -printf '%y %12s  %TY-%Tm-%Td %TH:%TM  %p\n' 2>/dev/null \
  | sort -k5 > "logs/manifest_$STAMP.txt"
echo "captured file inventory ($STAMP): $(wc -l < "logs/manifest_$STAMP.txt") entries"

# 4. housekeeping: keep logs/ lean. Retain the newest few ephemeral snapshots this
#    script writes (state_/manifest_) and cap the raw slurm-*.out it sweeps in; the
#    named job logs (logs/**/<name>.out) are the record and are left alone.
ls -t logs/state_*.txt    2>/dev/null | tail -n +6  | xargs -r rm -f
ls -t logs/manifest_*.txt 2>/dev/null | tail -n +6  | xargs -r rm -f
ls -t logs/slurm-*.out    2>/dev/null | tail -n +21 | xargs -r rm -f

# 5. reconcile BOTH ways, always. Commit only logs/ (a re-run job rewrites tracked
#    logs -> a dirty tree, which is what makes a plain `git pull` fail; committing
#    logs/ clears it). Untracked junk elsewhere is left for you to deal with, not
#    swept into a commit. Then ALWAYS rebase to bring Claude's scripts DOWN, push UP.
#    ==> Run THIS, never `git pull` by hand. It is the only sync command.
git add -A logs/
git commit -q -m "cluster sync $STAMP" 2>/dev/null || true
git pull --rebase -q origin main && echo "pulled latest (Claude's scripts are now here)"
git push -q origin main 2>/dev/null && echo "pushed logs -> git" || true
echo "synced both ways."
