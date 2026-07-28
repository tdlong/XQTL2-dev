# AGE_SY catalog v6 — filter decision

BAQ-on catalog (#23). Counted 44 BAMs at the loosest filter (min-dp 10, maxaf 0.07,
snpgap 0), then swept indel-distance x maxaf after the fact and scored concordance vs v3
(fraction of (SNP,sample) pairs whose ALT freq differs >2% / >5%).

## Result (genome-wide)
| indel>= | maxaf | SNPs kept | >2% off v3 | >5% off v3 |
|--------:|------:|----------:|-----------:|-----------:|
| 0  | 3% | 1,355,018 | 4.81% | 0.74% |
| 0  | 5% | 1,370,502 | 4.89% | 0.79% |
| 0  | 7% | 1,376,802 | 4.94% | 0.81% |
| 10 | 3% | 1,283,110 | 4.24% | 0.40% |
| 10 | 5% | 1,295,786 | 4.28% | 0.42% |
| 10 | 7% | 1,300,606 | 4.31% | 0.43% |
| 20 | 3% | 1,207,616 | 4.03% | 0.32% |
| 20 | 5% | 1,219,120 | 4.07% | 0.33% |
| 20 | 7% | 1,223,424 | 4.09% | 0.34% |
| 50 | 3% | 1,028,252 | 3.81% | 0.28% |
| 50 | 5% | 1,037,812 | 3.84% | 0.29% |
| 50 | 7% | 1,041,344 | 3.86% | 0.30% |

## Decision: min-dp 10, maxaf 3%, snpgap 20
- indel is the lever; 20 is the knee. 20->50 gains only 0.32%->0.28% (>5%) but loses 180K SNPs.
- maxaf is a coin flip (3->7% = +1.3% SNPs, concordance drifts 0.32->0.34% the wrong way); take 3%.
- The >2% column has a ~3.8% floor = finite-depth counting noise; the >5% column is the real signal.
- BAQ-on discovery already suppressed most near-indel junk, so even snpgap 0 is decent (0.74%).

## Open
- #24: catalog.tsv.gz still holds ~3K multiallelic positions until fixed; the workaround
  merge_dedup.R drops them. Real fix = drop multiallelic positions in the catalog.
