# CLAUDE.md — XQTL2-dev

## What this repo is

This is the **Long lab private project repo**. It contains project-specific
helpfiles, submission scripts, and analysis code for lab experiments.
The pipeline lives in a separate repo (XQTL2) accessible via the `pipeline/`
symlink in this directory.

## Rules — read before doing anything

0. **Read `pipeline/Slurm.md` before writing any SLURM job.**
   It defines partition memory limits, core allocation rules, and job
   templates for this cluster. Key fact: standard partition is 6 GB max
   per core — requesting more memory than one core allows implicitly
   pulls additional cores.

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

## Scientific tables

Tables in manuscripts and supplements follow this format. `temp_aging/tbl_fg3.gif`
is a worked example; `temp_aging/style_tables.py` applies it to pandoc's docx
output and `temp_aging/make_supplement.sh` is the build.

```
────────────────────────────────────────────────────────────   full width, heavy
                                  h²
                        ──────────────────────────             spans only the columns it modifies
                     n      2%      1%     0.5%
────────────────────────────────────────────────────────────   full width, under the column headers
  Outbred panel   1,000    0.35    0.04    0.004
  Inbred lines      200    0.02    0.00    0.000
  Eight-way RILs    750    0.70    0.16    0.022
  X-QTL         114,000    1.00    1.00    0.894
────────────────────────────────────────────────────────────   full width, heavy

Note. — ...
```

- Never a vertical rule. Never a double rule.
- Four rules only: heavy at top, under the column headers, heavy at the bottom,
  and a lighter one under any superheader spanning **only** the columns it
  modifies.
- Rules are graded in weight: heaviest top and bottom, lighter under the column
  headers, lightest for a superheader.
- The stub column takes no header; row labels sit in it.
- Numbers right-aligned, column headers and superheaders centred over them.
- A `Note. —` sits below the bottom rule.

Two traps when generating this with pandoc:

- pandoc's default `reference.docx` defines a `Table` style with
  `<w:tblStylePr w:type="firstRow">` that draws a full-width rule under row 1 —
  the superheader in a grouped table. Strip it, or cell borders are overridden.
- Never put math (`$h^2$`) in a table header. It becomes an equation, which
  renders blank in some viewers and cannot be centred or spanned as text. Use
  plain text, `h²`.
