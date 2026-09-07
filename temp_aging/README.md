# temp_aging — AGE_SY longevity X-QTL manuscript

**The Google Doc "Longevity" is canonical.** This folder feeds it; it does not
feed back. When the Doc and a file here disagree, the Doc is right.

## Root — the manuscript, and nothing else

| file | what it is |
|---|---|
| `METHODS.md` | computational Methods |
| `RESULTS.md` | Results |
| `DESIGN_POWER.md` | Supplementary Result X — power of four designs |
| `Methods.docx` | built from `METHODS.md`, for pasting into the Doc |
| `Supplement.docx` | built from the supplement sources, likewise |

Rebuild both with `bash temp_aging/build/make_docs.sh`. Close them in Word
first or the build writes under an open file.

## Folders

- `legends/` — the four figure legends (`Figure1`, `Figure2`, `FigureS1`, `FigureS2`)
- `analysis/` — every `.R` that produces a number or a figure, the figures
  themselves, and `TableS1.tsv`. Drivers: `reproduce.sh`, `make_figures.sh`
  (fetches from HPC3, then draws), `run_numbers.sh` (runs each analysis and
  stamps its output with input MD5s and the pipeline commit).
- `build/` — document machinery only: `make_docs.sh`, `make_supplement.sh`,
  `reference.docx`, `fix_math.py` (repairs pandoc's `\bar` accents),
  `style_tables.py` (booktabs rules and column widths), `tbl_fg3.gif` (the
  table style we are matching).
- `numbers/` — output of `run_numbers.sh`; every number in the prose is here
  with the provenance of the data that produced it.
- `superseded/` — dead. Nothing in the manuscript depends on it.

## Standing constraints

Never push without being asked. HPC3 is unreachable from the assistant side —
everything comes back through git.
