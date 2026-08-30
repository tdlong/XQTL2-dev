# Methods — sequencing, alignment, SNP calling

Draft. `[[ ]]` marks facts not recorded in this repo.

## Sequencing

Genomic DNA from each pool — one cage, treatment and sex — was prepared into a
Nextera library [[kit and catalogue number]] and sequenced paired-end
[[2 × 150 bp]] on an Illumina NovaSeq X Plus at the University of California,
Irvine Genomics Research and Technology Hub. Sixty libraries were sequenced: ten
replicate cages × three groups (unselected control, SY10 selected, SY20
selected) × two sexes.

## Alignment

Reads were aligned to the *Drosophila melanogaster* dm6 reference genome with
`bwa mem` (Li 2013), coordinate sorted and indexed. Sample alignments were not
deduplicated; the eight founder alignments were.

## SNP calling

SNPs were called on the five major chromosome arms, X, 2L, 2R, 3L and 3R. A
catalog of SNPs segregating among the eight founders of the synthetic population
was built once and used for every library. Founders were called jointly
with `bcftools`, and a site entered the catalog if it was biallelic; covered at
10× or more in every founder; at alternate allele frequency of 3% or less, or 97%
or more, in every founder; fixed for the reference allele in at least one founder
and for the alternate in at least one other; and at least 20 bp from the nearest
founder indel, a distance chosen from a threshold sweep on these data. Founder B5
was exempted from these criteria on chromosome 2L, where its coverage is
insufficient to evaluate them. Of 1,887,667 initial candidate positions,
1,207,436 entered the catalog [[confirm against
process/AGE_SY/Catalog/catalog.stats.txt]].

Reference and alternate read counts for each pooled library were then taken
directly from the pileup at each catalog position, giving the allele counts from
which haplotype frequencies were estimated.

## Coverage

Sequencing depth was summarised as each library's median reference plus
alternate read count across catalog SNPs. Across the sixty libraries,
autosomal depth averaged 140× per library, with per sex medians of 129× and
127× and a range of 49× to 400×; fifteen libraries fell below 75× and twenty
below 100×. On the X chromosome female libraries averaged 134× and male
libraries 72×, the latter 0.54 of the former, consistent with male
hemizygosity. No library fell below the pipeline's coverage or uniformity
thresholds.

## Haplotype estimation

The mapping population descends from eight inbred founders, so a pool's allele
frequency at any SNP is a mixture: the sum over founders of each founder's
frequency in the pool, weighted by its genotype at that SNP. The founder
genotypes are known, so with many SNPs in a local window the eight founder
frequencies can be recovered by least squares.

Founder frequencies were estimated at positions every 5 kb. Each position was
given a window of ±75 kb, widened in inverse proportion to the local
recombination rate so that windows span comparable genetic rather than physical
distance; recombination rate was taken from an eighth degree polynomial fitted
to each chromosome arm. Windows holding fewer than 50 catalog SNPs were dropped.
At each window, and for each library independently, the eight frequencies were
estimated by constrained least squares (`lsei`, R package `limSolve`),
minimising the squared deviation between observed and predicted allele frequency
subject to the frequencies summing to one and each being non-negative. The fit
also returns a covariance matrix for the estimates; this was averaged over the
ten windows either side of each position, ±50 kb, and carried into the tests.

Founders cannot always be told apart. Where two or more carry the same haplotype
across a window their individual frequencies are unidentifiable and only their
sum is determined. Founders were therefore clustered within each window by
hierarchical clustering of their genotypes, cut at a height of 2.5, and any
founder sharing a cluster with another was masked. Masked windows, together with
the few in which the fit failed to satisfy its own constraints, were filled by
linear interpolation between the means of the flanking resolved windows. The
filled series was then smoothed with a running mean spanning ±100 kb, that is
the twenty estimation windows either side of each position.

## Software

Alignment, SNP calling, haplotype estimation and the genome scans were performed
with the XQTL2 pipeline (github.com/tdlong/XQTL2), where tool versions and
parameters are documented.

---

**Still needed:** the Nextera kit and read length, and the number of flies per
pool and the DNA extraction method, neither recorded in this repo.

**Belongs in the Design section, not here:** that two further cages were run and
excluded as a separate source population. This section describes the ten
replicates the paper analyses.
