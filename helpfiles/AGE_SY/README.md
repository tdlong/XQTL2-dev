# helpfiles/AGE_SY

The analysis is **10 replicates from the November 2023 cage** — AGE_SY replicates
1–7 and 10–12. Replicates 8 and 9 came from the May 2023 cage and are excluded
(`helpfiles/AGE_2024/population_assignment.txt`). Those 10 split 5/5 into odd
(1,3,5,7,11) and even (2,4,6,10,12), which is where the replicate error comes
from.

## What the analysis uses

| file | what |
|---|---|
| `nov_only/<scan>.no89.txt` | the four 10-replicate designs — **the scans** |
| `nov_only/<scan>.no89.{odd,even}.txt` | the 5/5 halves — **the error term** |
| `AGE_SY_haplotype_parameters_size75k.R` | 75 kb windows, 5 kb step, 8 B founders |
| `AGE_SY.bams` | the 72 sample BAMs + 8 founders, for `call_samples.sh` |

Built by `scripts_oneoffs/AGE_SY/nov_only/make_designs.R` from the four
`AGE_SY<diet>_<sex>.test.txt` files, which are the full 12-replicate designs and
are kept as the source those are cut from.

**`Num` and `Proportion` in the design files are the authoritative fly counts.**
They are what the pipeline reads and what the heritability calculation uses. All
40 selected pools have both.

## barcodes/

Barcode-to-sample maps, one per sequencing run. History — used by `fq2bam` when
the reads were first aligned, not by anything downstream.

## retired/

Nothing here feeds the current analysis.

- `AGE_SY*_R1to6.test.txt` — first 6 replicates, for the AGE_2024 pilot comparison
- `AGE_SY*_{May,Nov}.test.txt` — split by source cage, for the same comparison
- `splithalf/` — odd/even over all 12 replicates; superseded by `nov_only/*.odd/even`
- `AGE_SY_haplotype_parameters.R` — 150 kb windows; the size75k file replaced it
- `info.AGE_SY.txt` — a fly-count summary with blank rows for R11/R12. **Do not
  use it for counts**; it was misread once already. The design files are the record.
- `temp.AGE_SY.txt`, `bam_list.v4*.txt` — referenced by nothing
