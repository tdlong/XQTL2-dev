# logs/ — cluster↔git observability bridge

Claude Code runs off-cluster and can read only what is in git. This directory is
the channel that gives it full visibility into cluster runs without anyone
pasting terminal output.

## How it works
- **Tracked, subdivided by project:** `logs/<project>/` (e.g. `logs/AGE_SY/`).
  `.gitignore` whitelists `logs/**`, overriding the repo-wide `slurm*` / `*.out`
  ignores, so log and state files here are committed.
- **Every SLURM script writes its `-o` here:** `#SBATCH -o logs/<project>/<name>_%j.out`.
- **Pipeline jobs** (whose `-o` we don't control) drop `slurm-<jobid>.out` in the
  repo root; `bin/cluster_sync.sh` sweeps those into `logs/<project>/`.

## The one command (run on the cluster after any activity)
```bash
bash bin/cluster_sync.sh <project> [process_dir]
#   sweeps slurm-*.out, captures sacct + squeue state, writes an output manifest,
#   then commits + pushes everything under logs/.
```
Then Claude runs `git pull` and reads `logs/<project>/` — logs, job state, and
what files the run produced.

## Housekeeping
These are just tracked text files. Delete old ones anytime (`rm logs/<project>/state_*.txt`)
and commit — no retention policy enforced.
