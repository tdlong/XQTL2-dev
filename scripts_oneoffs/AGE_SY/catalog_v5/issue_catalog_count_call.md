## Summary

`catalog_count.sh` takes each sample's per-SNP `FORMAT/AD` from **after**
`bcftools call -m -C alleles`. At a known catalog SNP where a sample carries a
low-but-real minor-allele fraction, `call -m` applies a **per-sample diploid
genotype model**, calls the site homozygous-REF, drops the ALT allele (`ALT="."`),
and the ALT read count **collapses to zero** — even though the ALT reads are
physically present and pass every mpileup quality filter. For pooled sequencing
this deletes real signal. **The counter should read `AD` straight from `mpileup`,
not from the genotype-collapsed record.**

The current recipe (`catalog_count.sh:53-55`):
```
bcftools mpileup -B -q20 -Q20 --max-depth 2000 -T "$cat" -a FORMAT/AD \
  | bcftools call -m -C alleles -T "$cat" -Ob > "$tmp"
```

## Evidence — 5 clean chrX SNPs, sample Con_R5_F (AGE_SY)

All 5 are **founder-clean** (segregating, all founders fixed, no intermediate,
`dist_indel` 41–497 bp so NOT near an indel). ALT reads counted by `mpileup`
after the full v4 quality filters (`-B -q20 -Q20`) vs after the `call` step:

| SNP | catalog REF,ALT | ALT reads from mpileup (`-B -q20 -Q20`) | ALT after `\| call -m -C alleles` |
|-----|-----------------|------------------------------------------|------------------------------------|
| chrX:6471332  | A,G | 33 | **0** (ALT→".") |
| chrX:12335553 | A,T | 26 | **0** |
| chrX:10242268 | A,C | 37 | **0** |
| chrX:10928272 | C,G | 35 | **0** |
| chrX:5332918  | C,A | 45 | 48 (kept — ~18% ALT clears the caller's threshold) |

The **only** difference between the two columns is the `| bcftools call -m -C alleles`
step. `-q20`, `-Q20`, and BAQ do not remove these reads (they survive to 33/26/37/35).
The genotype caller does. The one site that survived had a slightly higher ALT
fraction, i.e. it is a caller likelihood threshold, not a read filter.

## Reproduce (one site)
```
CHR=chrX; POS=6471332; REF=<dm6.fa>; CAT=<Catalog>/catalog.tsv.gz; BAM=<Con_R5_F.bam>
# mpileup AD — 33 ALT reads are there:
bcftools mpileup -B -q20 -Q20 --max-depth 2000 -r $CHR:$POS -T "$CAT" -a FORMAT/AD -f $REF $BAM \
  | bcftools query -f '[%AD]\n'          #  -> 212,33
# after the caller (the real catalog_count recipe) — ALT erased:
bcftools mpileup -B -q20 -Q20 --max-depth 2000 -r $CHR:$POS -T "$CAT" -a FORMAT/AD -f $REF $BAM \
  | bcftools call -m -C alleles -T "$CAT" | bcftools query -f '[%AD]\n'   #  -> 227
```

## Why v3 (the validated `bam2bcf2REFALT.sh`) does not have this
v3 runs `mpileup | call -mv` **jointly over all samples** (`-b bams`), so a
per-sample-rare allele is cohort-confident and each sample's AD is retained
(Con_R5_F reports 56/291 ALT at chrX:6471332). The catalog counter calls each
sample **alone**, so the diploid model sees ~13% ALT in isolation and discards it.

## Scale
For a single sample (Con_R5_F), chrX, mid-arm, coverage 50–300× on both callers:
3,537 SNPs disagree v3-vs-v4 by >10% ALT freq. Founder triage: 493 (14%) are
near-indel (a separate, legitimate filtering question); **3,044 (86%) are
founder-clean, and 2,961 of those are v4-loses-ALT** — this mechanism. Almost all
are `mpileup ALT present → post-call ALT = 0`.

## Fix
This is pooled sequencing at a **known** SNP — the task is to **count REF vs ALT
reads**, not to genotype. Remove the genotype model from the counting path: read
per-sample `AD` from `bcftools mpileup -a FORMAT/AD` at the catalog `-T` positions
**directly**, and do not route it through `bcftools call -m -C alleles`. The only
legitimate exclusions are bad *reads* — low base quality, base near the read end,
etc. — which `mpileup`'s own `-q`/`-Q`/read-position handling already applies. A
per-sample diploid genotype prior ("is this a variant?") does not apply to a pooled
sample and must not re-decide whether an already-catalogued SNP exists. (If `call`
is wanted only to order alleles to the catalog REF,ALT, take depths from the mpileup
record, not the GT-collapsed post-call record.)
