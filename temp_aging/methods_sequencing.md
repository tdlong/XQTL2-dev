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

**Coverage.** Depth was taken as each library's median REF+ALT count over
catalog SNPs. Across the 60 libraries the scans use, autosomal depth averaged
140× per library (per-sex medians 129× and 127×; range 49–400×). On chrX female
libraries averaged 134× and male libraries 72×, the latter 0.54 of the former, as
hemizygosity predicts; X and autosomes are therefore not pooled into a single
figure. No library was flagged low-coverage or patchy by the pipeline's QC.

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

Coverage is sourced: `temp_aging/coverage.R` -> `numbers/coverage.txt`, reading
`process/AGE_SY/Calls/refalt_qc.txt` from the pipeline's own `refalt_qc.R`. It
restricts to the 60 libraries the scans use, and reports chrX apart from the
autosomes because males are hemizygous. Both were run on HPC3 on 2026-08-29.

An earlier draft of this file quoted depth bins from
`logs/AGE_SY/compare_summary.out` instead. Those are real but describe something
else — per (SNP, sample) pooled over all 72 libraries, at SNPs shared between the
old and new callers, with male chrX folded into the autosomes. **Do not reinstate
them.** They put 76% of observations in 50–199×; the actual per-library medians
average 140× on the autosomes, so the old figures understated depth as well as
mixing the wrong samples.

**Needs Tony:** Nextera kit and read length; number of runs; flies per pool and
DNA extraction, which are nowhere in this repo; and whether to report 72
libraries (all called) or 60 (the ten replicates scanned) — the catalog and
counts were built on all 72, the scans drop replicates 8 and 9.
