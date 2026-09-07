# temp_aging — what is canonical

**The Google Doc "Longevity" is the manuscript.** It is the live copy of the
title, intro, Results, figure legends 1–3 and the wet-lab Methods. Prose here
feeds it; it does not feed back. When the two disagree, the Doc wins.

## What this directory is for

Analysis, numbers, figures, and the text that has to be generated from them.

| | |
|---|---|
| `RESULTS.md`, `METHODS.md` | source for the computational text. Superseded by the Doc once pasted. |
| `DESIGN_POWER.md` | Supplementary Result X, the power comparison |
| `Figure*_legend.md` | one legend per figure, including `FigureS1_legend.md` and `FigureS2_legend.md` |
| `TableS1.tsv` | Supplementary Table 1, from `make_TableS1.R` |
| `Figure*_plot.png` | figures, each from its `make_*.R` |

## Build

```
bash temp_aging/make_docs.sh
```

Produces exactly two files, both for uploading to Drive and opening with Google
Docs:

| | |
|---|---|
| `Methods.docx` | the computational Methods, whole |
| `Supplement.docx` | Table S1, the S1 and S2 legends, Supplementary Result X |

Supporting, not run directly: `make_supplement.sh`, `style_tables.py`
(scientific table format, see CLAUDE.md), `fix_math.py` (repairs pandoc's bar
accents), `reference.docx` (1 in margins, no conflicting table style),
`make_power_table_html.sh` (the power table as HTML, when a paste works better).

`.docx` output is gitignored — regenerate it, do not edit it.

## superseded/

Retired drafts and bookkeeping from earlier phases: `MANIFEST.md`,
`NUMBERS_TO_UPDATE.md`, `PIPELINE.md`, `METHODS_yours.md`, `peak_table.txt`,
`make_maintext.sh` and older `METHODS.md` / results drafts. Kept for history,
not read by anything.
