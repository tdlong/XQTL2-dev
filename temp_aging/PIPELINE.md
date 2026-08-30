# PIPELINE — what produces what, and where it runs

Read this before rerunning anything. `MANIFEST.md` says what each file *is*;
this says what *makes* it and what *breaks* when it changes. Every path is
relative to the repo root.

`HPC3` = must run on the cluster, `local` = runs on the laptop from files
fetched down. Nothing in the HPC3 column can run here: `process/` is not on
this machine.

## Two primary products, everything else derived

A scan writes exactly two things that are not a function of something else:

    Scans/<scan>/<scan>.scan.txt                 Wald and h2 per window
    Scans/<scan>/<scan>.meansBySample.<chr>.txt  founder frequencies per pool

Every table, number and figure below is a function of those two. So the question
that matters is not "what is each file" but "given a new scan, how do I rebuild
everything". The answer is one command:

    sbatch temp_aging/reproduce.sh      # on HPC3, from the repo root

It gathers, partitions, subsets the zoom peaks and reruns the numbers, in order,
reporting which step failed rather than dying silently partway. It resubmits
nothing. `run_scans.sh` submits it as the dependent job after the twelve scans,
so there is one definition of "derived", used both ways.

## Getting to the primary products

| # | Stage | Runs | Command | Writes |
|---|---|---|---|---|
| 1 | align | HPC3 | `pipeline/scripts/fq2bam.sh` | `data/bam/AGE_SY/*.bam` |
| 2 | catalog | HPC3 | `pipeline/scripts/build_catalog.sh` | `process/AGE_SY/Catalog/` |
| 3 | count | HPC3 | `pipeline/scripts/call_samples.sh` | `Calls/RefAlt.<chr>.txt` |
| 4 | sample QC | HPC3 | `pipeline/scripts/refalt_qc.R` | `Calls/refalt_qc.txt` |
| 5 | haplotypes | HPC3 | `pipeline/scripts/run_haps.sh` | `Haps/R.haps.<chr>.out.rds` |
| 6 | **12 scans** | HPC3 | `nov_only/run_scans.sh` | **the two primary products** |

## Derived, all of it by reproduce.sh

| Product | From | Consumed by |
|---|---|---|
| `AGE_SY_4scan_no89.txt.gz` | the 4 scan tables | Fig 1a, 1b; Fig 2 Wald; Fig rr; all five numbers scripts |
| `AGE_SY_splithalf_H2_no89.txt.gz` | the 8 half-scan tables | `significant_regions.R`; input to the partition |
| `H2_varcomp_by_window_no89.txt.gz` | the split-half table | **Fig 1c only** |
| `AGE_SY_zoom_means.txt.gz` | the means files | Fig 2; `chr3L_peak.R` |
| `numbers/*.txt`, `peak_table.txt` | the three above | the prose |

## Outside reproduce.sh

`AGE_SY_4snpscan_no89.txt.gz` (Fig 3) is the exception. It is imputed from the
**smoothed haplotypes**, so a new scan makes it stale, but it is a separate
submission and nothing enforces the order:

    bash scripts_oneoffs/AGE_SY/nov_only/run_snp_scans.sh

It takes no `--sex` of its own and must not be given one. The chrX dosage is
already in the `Num` it reads: `smooth_haps.R:200` scales `Num` before writing
the rds, and `snp_scan.R` pulls `Num` from that rds, never from the design file.
Applying `--sex` again would square the factor. Which sex an rds was built under
is recorded -- `readRDS(f)$sex` -- so a stale SNP scan is detectable rather than
inferred (XQTL2 #39).

Figures are the other exception: they run on the laptop, from files fetched down.

## The trap inside the derived chain

`varcomp_H2.R` takes input and output as arguments and its defaults are the
*12-replicate* paths. Run bare it rebuilds the wrong file and Figure 1c goes on
reading a stale no89 partition, with no error anywhere. `reproduce.sh` always
passes both. Do not call it by hand without them.

## Figures

One command on the laptop, on VPN, from the repo root:

    bash temp_aging/make_figures.sh

It rsyncs the six `process/` files the figure scripts read -- one connection,
not six scp lines -- then draws all four figures and lists what it wrote.
`--no-fetch` to draw from what is already local; a trailing `1 2 3 rr` to draw
only some.

`numbers/*.txt` and `peak_table.txt` are NOT fetched: `reproduce.sh` writes them
on HPC3 and they come back through git.

## What --sex changes

`--sex` (XQTL2 #38) scales pool size on **chrX only**; autosomes get 1.0
whatever it says. So after a rerun with it, autosomal values must be identical
to the byte. That is the check on any rerun: if an autosomal number moves,
something other than the sex flag changed.

It enters at step 6, so everything from 6 downward is stale until rerun —
including the SNP scan, which is not chained. The SNP scan picks the correction
up through the rds, so rerunning it is sufficient; it needs no flag of its own.
