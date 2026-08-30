# PIPELINE — what produces what, and where it runs

Read this before rerunning anything. `MANIFEST.md` says what each file *is*;
this says what *makes* it and what *breaks* when it changes. Every path is
relative to the repo root.

`HPC3` = must run on the cluster, `local` = runs on the laptop from files
fetched down. Nothing in the HPC3 column can run here: `process/` is not on
this machine.

## The chain

| # | Stage | Runs | Command | Writes |
|---|---|---|---|---|
| 1 | align | HPC3 | `pipeline/scripts/fq2bam.sh` | `data/bam/AGE_SY/*.bam` |
| 2 | catalog | HPC3 | `pipeline/scripts/build_catalog.sh` | `process/AGE_SY/Catalog/` |
| 3 | count | HPC3 | `pipeline/scripts/call_samples.sh` | `process/AGE_SY/Calls/RefAlt.<chr>.txt` |
| 4 | sample QC | HPC3 | `pipeline/scripts/refalt_qc.R` | `process/AGE_SY/Calls/refalt_qc.txt` |
| 5 | haplotypes | HPC3 | `pipeline/scripts/run_haps.sh` | `process/AGE_SY/Haps/R.haps.<chr>.out.rds` |
| 6 | **12 scans** | HPC3 | `scripts_oneoffs/AGE_SY/nov_only/run_scans.sh` | `process/{AGE_SY,AGE_SY_splithalf}/Scans/<scan>/` |
| 7 | gather | HPC3 | chained on 6 | `AGE_SY_4scan_no89.txt.gz`, `AGE_SY_splithalf_H2_no89.txt.gz` |
| 8 | partition | HPC3 | chained on 7 | `AGE_SY_splithalf/H2_varcomp_by_window_no89.txt.gz` |
| 9 | zoom means | HPC3 | chained on 6 | `AGE_SY/AGE_SY_zoom_means.txt.gz` |
| 10 | SNP scans | HPC3 | `scripts_oneoffs/AGE_SY/nov_only/run_snp_scans.sh` | `AGE_SY/AGE_SY_4snpscan_no89.txt.gz` |
| 11 | numbers | HPC3 | chained on 7 and 9 | `temp_aging/numbers/*.txt` |
| 12 | figures | local | `temp_aging/make_figure*.R` | `figures/*.png` |

Steps 6-9 and 11 are one command: `run_scans.sh` submits them all with SLURM
dependencies and returns. Step 10 is a **separate submission** — it reads the
smoothed haplotypes step 6 writes, so it must be run after, and it is not
chained. Step 12 needs the files fetched down first.

## What each product feeds

| Product | Consumed by |
|---|---|
| `AGE_SY_4scan_no89.txt.gz` | Fig 1a, 1b; Fig 2 Wald panels; Fig rr; **all five** numbers scripts |
| `AGE_SY_splithalf_H2_no89.txt.gz` | `significant_regions.R`; input to step 8 |
| `H2_varcomp_by_window_no89.txt.gz` | **Fig 1c only** |
| `AGE_SY_zoom_means.txt.gz` | Fig 2; `chr3L_peak.R` |
| `AGE_SY_4snpscan_no89.txt.gz` | **Fig 3 only** |
| `Calls/refalt_qc.txt` | `coverage.R` |

## The two that are easy to miss

**Step 8 is not gather.** `varcomp_H2.R` takes two arguments, input and output,
and its defaults are the *12-replicate* paths. The no89 partition only exists
because someone passed both:

    Rscript scripts_oneoffs/AGE_SY/splithalf/varcomp_H2.R \
        process/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz \
        process/AGE_SY_splithalf/H2_varcomp_by_window_no89.txt.gz

Run without arguments it silently rebuilds the 12-replicate file instead, and
Figure 1c goes on reading a stale no89 file with no error anywhere.

**Step 10 depends on step 6 but is not chained to it.** The SNP scan imputes
per-SNP frequencies *from the smoothed haplotypes*, so rerunning the haplotype
scans makes the SNP scan stale. Nothing enforces this.

## Fetching to the laptop

Every figure input, from the repo root, on VPN:

    R=tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2-dev/process
    mkdir -p process/AGE_SY/Calls process/AGE_SY_splithalf
    scp $R/AGE_SY/AGE_SY_4scan_no89.txt.gz                   process/AGE_SY/
    scp $R/AGE_SY/AGE_SY_zoom_means.txt.gz                   process/AGE_SY/
    scp $R/AGE_SY/AGE_SY_4snpscan_no89.txt.gz                process/AGE_SY/
    scp $R/AGE_SY/Calls/refalt_qc.txt                        process/AGE_SY/Calls/
    scp $R/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz  process/AGE_SY_splithalf/
    scp $R/AGE_SY_splithalf/H2_varcomp_by_window_no89.txt.gz process/AGE_SY_splithalf/

`numbers/*.txt` do NOT need fetching — `run_numbers.sh` runs on HPC3 and they
come back through git.

## What --sex changes

`--sex` (XQTL2 #38) scales pool size on **chrX only**; autosomes get 1.0
whatever it says. So after a rerun with it, autosomal values must be identical
to the byte. That is the check on any rerun: if an autosomal number moves,
something other than the sex flag changed.

It enters at step 6, so everything from 6 downward is stale until rerun —
including step 10, which is not chained.

`snp_scan.R` has **no** X dosage handling at all, not even the old fixed 0.75,
so Figure 3 on chrX is wrong independently of this fix. Rerunning step 10 does
not correct it.
