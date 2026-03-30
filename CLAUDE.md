# CLAUDE.md — XQTL2-dev

## What this repo is

This is the **Long lab private project repo**. It contains project-specific
helpfiles, submission scripts, and analysis code for lab experiments.
The pipeline lives in a separate repo (XQTL2) accessible via the `pipeline/`
symlink in this directory.

## Rules — read before doing anything

1. **Read the pipeline README before writing any script.**
   The pipeline is at `pipeline/` (symlink to the XQTL2 install).
   Read `pipeline/README.md` to understand what scripts already exist.

2. **Never rewrite pipeline scripts.** If a script exists in `pipeline/scripts/`,
   call it — do not copy or rewrite it here. Pipeline improvements belong in
   XQTL2, not here.

3. **All pipeline calls use `pipeline/scripts/` paths.** For example:
   ```bash
   bash pipeline/scripts/run_scan.sh --design helpfiles/ZINC2/design.txt ...
   ```

4. **Pipeline reference files are at `pipeline/helpfiles/`.** For example:
   `pipeline/helpfiles/het_bounds.txt`, `pipeline/helpfiles/founder.bams.txt`.
   Do not duplicate these here.

5. **Project data lives in `data/` and `process/`** — both gitignored.
   Never commit data files, BAMs, or pipeline output.

6. **This repo has one remote: `origin` → XQTL2-dev.git.**
   Never add XQTL2.git as a remote. Never push to the public pipeline repo.

7. **scripts_oneoffs/<project>/** contains submission scripts for that project.
   These call `pipeline/scripts/` and reference `helpfiles/<project>/`.
   No shared utilities — if something is used by more than one project it
   belongs in the pipeline (XQTL2), not here.

8. **helpfiles/<project>/** contains project configs: barcodes, design files,
   haplotype parameters. One subdirectory per project.

## Directory structure

```
XQTL2-dev/
├── pipeline -> ../XQTL2        ← symlink to pipeline (not tracked in git)
├── helpfiles/
│   └── <project>/              ← barcodes, design, hap params
├── scripts_oneoffs/
│   └── <project>/              ← submission scripts calling pipeline/scripts/
├── data/                       ← gitignored (raw reads, BAMs, founders)
├── process/                    ← gitignored (pipeline outputs)
├── figures/                    ← gitignored (publication figures)
└── logs/                       ← gitignored (SLURM logs)
```

## Updating the pipeline

When XQTL2 releases updates:
```bash
cd ../XQTL2
git pull origin main
```
Your scripts automatically use the new version via the symlink.
