# temp_aging — what is here, and what it rests on

> **Rebuilding after a new scan: see `PIPELINE.md`, not this file.** A scan
> produces two primary things -- the Wald/h2 tables and the means -- and
> everything else is derived from them by `bash temp_aging/reproduce.sh`. This
> file is an inventory of what each artifact *is*; it is the wrong place to look
> for what *makes* it.

AGE_SY lifespan XQTL: the Results prose, the figures, and the scripts that
produce every number in both. Written so nothing here has to be
reverse-engineered.

## The dataset

**Ten replicates.** Replicates 8 and 9 are dropped throughout — they are the May
2023 source cage, the rest are November 2023
(`helpfiles/AGE_2024/population_assignment.txt`). Dropping the pair keeps the
odd/even split balanced at 5 and 5. Every file used here has `no89` in its name;
the 12-replicate files still exist on disk and nothing current reads them.

    process/AGE_SY/AGE_SY_4scan_no89.txt.gz              the four scans
    process/AGE_SY/AGE_SY_4snpscan_no89.txt.gz           the four SNP scans
    process/AGE_SY/AGE_SY_zoom_means.txt.gz              founder freqs at 7 peaks
    process/AGE_SY_splithalf/AGE_SY_splithalf_H2_no89.txt.gz   odd/even halves

Produced on HPC3 by `scripts_oneoffs/AGE_SY/nov_only/` — `make_designs.R`,
`run_scans.sh`, `run_snp_scans.sh` (both chain their own gather).

## The unit

**268 one cM tiles.** Each arm is cut at 0–1, 1–2, 2–3 cM from its own start. A
tile is kept if any part of it is euchromatic by
`pipeline/helpfiles/het_bounds.txt`; 270 tiles, 2 wholly heterochromatic dropped,
5 straddling kept. Every genome-wide proportion in the prose is a proportion of
these 268.

**Thresholds are 7.5 and 15**, not 5 and 10. A tile counts as significant if some
window in it reaches 7.5 in some treatment: 98 of 268, 37%. Above 15: 12%.

**Support intervals have nothing to do with tiles.** They are read physically off
the scan — walk left and right from a peak while within 2 of it — and then
converted to cM.

## The prose

Tony's current text, stored clean with no annotations in the prose itself.

**Every number in them predates the h2 and covariance fixes of 30 August
(XQTL2 #40) and none should be quoted until the rerun lands.** The prose is
current; the figures in it are not. `NUMBERS_TO_UPDATE.md` lists each quoted
number, whether it survives, and which `numbers/` file replaces it -- including
the two places where the argument rather than the value has to change.

    RESULTS.md              all four paragraphs, in order

Paragraph 5 onward is not written.

    methods_sequencing.md   sequencing, alignment, SNP calling, coverage. NOT
                            signed off. Fills what METHODS.md never covered.
                            Written about the ten replicates only -- 60 libraries.
                            Runs and re-sequencing are deliberately absent: they
                            are logistics, and the depth numbers say what we got.
                            Needs the Nextera kit, read length, flies per pool and
                            extraction method -- none in this repo.

## The scripts, and the numbers they own

`run_numbers.sh` runs all five into `numbers/<script>.txt`, each with a header
naming every input it read, that input's size, mtime and MD5, and the git commit
of the script. **Any number quoted in the prose should be findable in one of
those files.** If it is not, it is not sourced.

    h2_threshold.R        the floor (tiles below Wald 2) and the fitted h2 at
                          7.5 and 15. Not a conversion -- h2 and Wald are
                          correlated but not monotone; these are averages
                          through wide scatter.
    scan_resolution.R     the tiling itself: 268 tiles, 264 cM, 122 Mb, and the
                          proportion above each threshold.
    peak_table.R          Table S1. Peaks top-down on the max across traits, 5 cM
                          exclusion, split at 15. Writes TableS1.txt, the supplementary peak table.
    significant_regions.R sex and diet over the 98 significant tiles. Each tile
                          is taken at ONE window, its peak across the four
                          treatments, so the eight split-half values describe the
                          same position.
    coverage.R            depth at catalog SNPs over the 60 libraries the scans
                          use. Reads Calls/refalt_qc.txt, which the pipeline's
                          own refalt_qc.R writes -- RUN BOTH ON HPC3, process/ is
                          there. chrX is kept apart from the autosomes: males are
                          hemizygous, so a single genome-wide depth would hide it.
    chr3L_peak.R          the 3L locus: the scan there, the founder shifts, and
                          the founder frequencies in the unselected controls
                          across the ten cages in order.

### Two things that are easy to get wrong

**The error term is the split halves, not a floor.** Five odd and five even
replicates give two independent h2 estimates of the same tile; half their squared
difference is a pure error term, and it is subtracted from every component of the
partition. The null-window floor in `h2_threshold.R` exists only because a single
h2 estimate has nowhere else to get an error term.

**Coverage bins in the old logs are not the coverage.**
`logs/AGE_SY/compare_summary.out` carries depth bins that put 76% of
observations in 50-199x. They are per (SNP, sample) pooled over all 72 libraries,
at SNPs shared between the old and new callers, with male chrX folded into the
autosomes -- a caller-comparison by-product. The real per-library medians average
140x on the autosomes. An early draft of methods_sequencing.md quoted the bins;
do not reinstate them.

**The partition's error term survives the h2 fix; the per-pool bias correction
does not.** Two different corrections were being applied and only one is
superseded. `h2_rep` (XQTL2 #40) removes the BIAS -- the inflation from squaring
noisy frequency estimates -- which is what subtracting `Cutl_H2_bias` per pool
used to do, badly. That path is now redundant and must not also be applied. The
split-half term removes the VARIANCE of the h2 estimate: subtracting `var(d)/n`
inside a window says nothing about how much that window's h2 would move on a
rerun, which is what the odd-even difference measures. It stays.

Expect the fractions to move anyway, for a reason that is not biology. The old
bias was ~0.5-1.2 and common to both halves, so it never entered the odd-even
difference -- it went into the SHARED component, 85.3%. Removing it from every
cell should cut shared and mechanically raise the sex (15.3%, then 13.3%) and
diet (1.6%) fractions. Do not compare those numbers across the change as though
they were like for like.

Open: whether to run the partition on `h2_rep` per half -- 5 replicates, so the
`var/n` correction is noisier -- or restructure it. This is ours to decide;
`varcomp_H2.R` is in `scripts_oneoffs/`, not the pipeline.

**The bias is not common to the four treatments.** `Cutl_H2_bias` is computed per
pool per window from that pool's own lsei reconstruction error and the
multinomial sampling of its own flies (XQTL2 #34). It runs 0.73 in SY10 females
to 0.83 in SY20 males, tracking how many flies were selected and how well each
pool reconstructs, so a sex contrast on raw h2 contains a sex difference in bias.
It is subtracted per pool per window before anything is contrasted. Leaving it in
halves every fraction in paragraph 3 (15.3% becomes 6.5%).

## The figures

PNGs are gitignored repo-wide (`*.png`), so they are not in git and have to be
regenerated. All read from `process/`; `make_figure3.R` scps its own input if it
is absent.

    make_figure1.R      3 panels. X_UNIT=cM gives the genetic axis as
                        Figure1_cM_plot.png, concatenated by linkage group.
    make_figure2.R      7 zoom cells plus a legend/control cell. DROP_REP c(8,9).
    make_figure3.R      SNP scan, females and males, both diets overlaid, y split
                        at 40, alpha a function of Wald.
    make_figure_rr.R    cM/Mb with mean Wald overlaid; why the base of 3L and 3R
                        is broad in Mb and not in cM.
    make_zoom_means.R   RUN ON HPC3. Subsets the 257 MB meansBySample files down
                        to 1.2 Mb around each of seven peaks.

    Figure1_legend.md, Figure2_legend.md   Tony's text, verbatim. Do not append.

## Superseded — do not quote

Moved to `superseded/`, with a README saying why. Nothing there is read by any
script, and no number in it is current.

`archive/` holds the AGE_2024 pilot figures and superseded scripts. It is
gitignored, so it exists only on this machine.

## Open

- Paragraph 5 onward not drafted.
- `METHODS.md` still describes 12 cages, 75 kb tiles and threshold 5, and its
  SNP calling describes the joint-QUAL caller that the founder-catalog caller
  replaced. `methods_sequencing.md` covers the calling; the rest is unrewritten.
- Whether Table S1 becomes a supplement is undecided.
- The power claim in paragraph 2 — that a 1–2% locus would have gone undetected
  in several hundred RILs — is asserted, not calculated. No design was named to
  compute it against.
- The chr3L support interval is **chr3L:9,280,000–9,335,000**, 55 kb, 0.219 cM
  (28.091–28.310). What genes are in it has not been checked. An earlier draft
  named the PGRP complex; that name has never been verified against an
  annotation and is not in any current file.
