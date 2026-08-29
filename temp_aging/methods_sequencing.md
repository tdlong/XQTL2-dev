# Methods — libraries, alignment, SNP calling, coverage

> DRAFT. Fills the gap in the superseded `METHODS.md`. `[[ ]]` items are not in
> this repo and need Tony. Detail beyond this is in the XQTL2 manual; the point
> here is what was done and what was chosen. Coverage caveat at the bottom.

**Libraries and sequencing.** Each pooled sample — one cage × treatment × sex —
was prepared as a Nextera library [[kit + catalogue no.]] and sequenced
paired-end [[2 × 150 bp?]] at the UCI Genomics Research and Technology Hub on an
Illumina NovaSeq X Plus. Seventy-two libraries (12 cages × 3 groups × 2 sexes)
were sequenced over [[N]] runs; replicates 8–11 were sequenced twice and merged
at the BAM level.

**Alignment.** `bwa mem` (Li 2013) to the *D. melanogaster* dm6 reference, five
major arms only (chrX, 2L, 2R, 3L, 3R). Pooled libraries are **not**
deduplicated: at 50–200× from Nextera, reads sharing a start site are real
molecules, and marking them as duplicates would delete them from the allele
frequencies the design measures. The eight founder BAMs are deduplicated, being
single inbred genotypes.

**SNP calling.** Pooled samples are never genotyped. A SNP catalog is built once
from the eight founders, and each sample is then counted at exactly those
positions, REF and ALT read depths taken directly from the pileup. Running a
diploid genotype caller on a pool would call a rare minor allele homozygous
reference and erase its reads, so none is used.

A founder SNP entered the catalog if it was covered at ≥ 10× in every founder,
sat at ALT frequency ≤ 3% or ≥ 97% in each of them, segregated (at least one
founder fixed for each allele), and lay ≥ 20 bp from the nearest founder indel —
that last threshold set by an indel-distance sweep on these data. Founder B5's
coverage collapses on chr2L, so B5 was exempted from these rules on that arm
alone. This gave **1,207,436 catalog SNPs** from 1,887,667 candidate positions,
the largest single loss (294k) being proximity to an indel.

**Coverage.** [[NOT YET SOURCED — run `temp_aging/coverage.R`. Wanted: mean
median-depth per library at catalog SNPs over the 60 kept libraries, autosomes
and chrX reported apart.]]

**Code.** Alignment, calling and everything downstream used the XQTL2 pipeline
(github.com/tdlong/XQTL2), which documents the tool versions and parameters.

---

## Provenance

Read out of the pipeline and logs, and safe to quote: tool versions and every
flag (`pipeline/scripts/fq2bam.sh`, `catalog_build.sh`, `catalog_count.sh`,
`catalog_filter.sh`); dm6 (`pipeline/ref/dm6.fa`); the no-deduplication rationale
(`pipeline/scripts/bam_qc.sh`, `pipeline/README.md`); the catalog tally and
thresholds (`logs/AGE_SY/clean_v6_filter.out`); B5's chr2L depth collapse
(`logs/AGE_SY/compare_v3_v4_54375203.out`); 72 sample + 8 founder BAMs
(`helpfiles/AGE_SY/AGE_SY.bams`).

**Coverage is not sourced yet.** An earlier draft of this file quoted depth bins
from `logs/AGE_SY/compare_summary.out`. Those are real, but they are per (SNP,
sample) pooled over all 72 libraries at SNPs *shared between the old and new
callers* — a caller-comparison by-product. Wrong sample set (includes replicates
8 and 9), wrong SNP set (blind to the ~26k chr2L SNPs only the current caller
has), and no separation of male chrX. They have been removed; do not reinstate
them.

The right number is the median depth per library at catalog SNPs, over the 60
libraries the scans actually use, with chrX apart from the autosomes because
males are hemizygous. `pipeline/scripts/refalt_qc.R` computes the per-sample
depths; `temp_aging/coverage.R` restricts them to the 60 and prints the summary.
Neither has been run against the current `process/AGE_SY` — refalt_qc.R must run
on HPC3 and the table be fetched back. The commands are in the header of
`coverage.R`, which is registered in `run_numbers.sh`, so once the input is in
place the number lands in `numbers/coverage.txt` like every other.

**Needs Tony:** Nextera kit and read length; number of runs; flies per pool and
DNA extraction, which are nowhere in this repo; and whether to report 72
libraries (all called) or 60 (the ten replicates scanned) — the catalog and
counts were built on all 72, the scans drop replicates 8 and 9.
