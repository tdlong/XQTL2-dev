#!/usr/bin/env bash
# cluster_sync.sh — the ONE command to run after any cluster activity.
#
# Claude Code runs off-cluster and can see only what is in git. This bridges the
# gap: it captures everything Claude needs to understand a run — the job logs,
# the SLURM job state, and a manifest of the output dir — into the tracked
# logs/<project>/ tree, then commits and pushes. Claude then `git pull`s and
# reads it. No pasting terminal output, ever.
#
# General (any project), not per-experiment:
#   bash bin/cluster_sync.sh <project> [process_dir]
#   e.g.  bash bin/cluster_sync.sh AGE_SY process/AGE_SY_v4
#
# Run from the repo root on the cluster. Old logs under logs/<project>/ can be
# deleted anytime (they are just tracked text).
set -uo pipefail
PROJ=${1:?usage: cluster_sync.sh <project> [process_dir]}
PDIR=${2:-}
LOGDIR="logs/$PROJ"
mkdir -p "$LOGDIR"
STAMP=$(date '+%Y%m%d_%H%M%S')

# 1. sweep the pipeline's default logs (slurm-<jobid>.out, dropped in CWD by jobs
#    we don't control the -o for) into the tracked tree.
shopt -s nullglob
n=0; for f in slurm-*.out; do mv "$f" "$LOGDIR/"; n=$((n+1)); done
[[ $n -gt 0 ]] && echo "swept $n slurm-*.out -> $LOGDIR"

# 2. capture SLURM job state: recent history + current queue.
{
  echo "### sacct (since yesterday)  $STAMP"
  sacct -X --starttime "$(date -d 'yesterday' +%F 2>/dev/null || date -v-1d +%F)" \
     --format=JobID,JobName%20,State,ExitCode,Elapsed,Start,End 2>/dev/null
  echo; echo "### squeue -u $USER"
  squeue -u "$USER" 2>/dev/null
} > "$LOGDIR/state_$STAMP.txt"
echo "wrote $LOGDIR/state_$STAMP.txt"

# 3. manifest of the output dir (names/sizes, not the bulk data) so Claude sees
#    what was produced without the files themselves (they stay in gitignored process/).
if [[ -n "$PDIR" && -d "$PDIR" ]]; then
  { echo "### manifest $PDIR  $STAMP"; ls -laR "$PDIR" | head -500; } > "$LOGDIR/manifest_$STAMP.txt"
  echo "wrote $LOGDIR/manifest_$STAMP.txt"
fi

# 4. one commit + push.
git add logs/
if git diff --cached --quiet; then echo "nothing new to sync"; exit 0; fi
git commit -q -m "cluster sync: $PROJ $STAMP"
git push -q && echo "synced -> git.  Claude: git pull && read $LOGDIR/"
