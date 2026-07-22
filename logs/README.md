# logs/ — cluster↔git observability bridge

Claude Code runs off-cluster and can read only what is in git. This directory is
how it sees cluster runs without anyone pasting terminal output.

## The one command (always identical, no arguments)
After any cluster activity, on the cluster, from the repo root:
```bash
bash scripts/cluster_sync.sh
```
It sweeps pipeline `slurm-*.out` logs, captures `sacct` + `squeue` state and a
manifest of everything under `process/`, writes them here, and commits + pushes.
Claude then `git pull`s and reads `logs/`.

## Details
- `logs/` is TRACKED; `.gitignore` whitelists `logs/**`, overriding the repo-wide
  `slurm*` / `*.out` ignores.
- SLURM scripts point their `-o` under `logs/` (e.g. `logs/AGE_SY/<name>_%j.out`),
  so those logs are already here; `cluster_sync.sh` adds the pipeline logs it
  can't set `-o` for, plus job state and the output manifest.
- Just tracked text. `cluster_sync.sh` auto-prunes its own `state_*`/`manifest_*`
  snapshots after 14 days; delete old job logs by hand whenever, then commit.
