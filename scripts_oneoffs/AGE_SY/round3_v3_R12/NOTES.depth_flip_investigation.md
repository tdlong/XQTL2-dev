# Depth-flag SNP-flip investigation — lab notebook

**Question:** in the RefAlt caller, changing `bcftools mpileup -d` (max reads
PER SAMPLE) from 1000 to 50000 keeps some SNPs and drops others. Why? And why do
the tiled (`-r` 5 Mb) and non-tiled (`-t` whole-chr) callers disagree at many
SNPs at the *same* `-d 1000`?

Keep this file updated every step so we do NOT re-run the hour-plus jobs.

---

## Pipeline facts (do not re-derive)

- Caller (`bam2bcf2REFALT.tiled.sh:54`, `bam2bcf2REFALT.sh:20`):
  `bcftools mpileup -I -d <D> -r <region> -a FORMAT/AD,FORMAT/DP -f ref -b bams | bcftools call -mv`
  Hardcoded `-d 1000`.
- The ONLY downstream filter (`bam2bcf2REFALT.tiled.sh:63`):
  `bcftools view -m2 -M2 -v snps -i 'QUAL>59'`  → so "dropped" = fails biallelic, SNP, or QUAL>59.
- `QUAL` is bcftools' POOLED site quality over all 80 samples. It is NOT a function of the per-sample AD counts.
- Modules: load **bcftools/1.21 ONLY**. Do NOT `module load samtools/1.10` in the
  same shell — it pulls htslib/1.10 and breaks bcftools (`undefined symbol
  bcf_has_variant_types, HTSLIB_1.16`). Load samtools in an isolated `( module
  purge; module load samtools/1.10; ... )` subshell.
- Data: `helpfiles/AGE_SY/AGE_SY.bams` = **80** samples. Ref `pipeline/ref/dm6.fa`.

## Cluster result locations (already computed — reuse, don't recompute)

- `process/prove_flips/call.d1000.vcf.gz`, `call.d50000.vcf.gz` — prove_depth_flips, 500 kb window
- `process/prove_flips/dc.d1000.vcf.gz`, `dc.d50000.vcf.gz` — depth_cap_check, ±10 kb windows
- `process/prove_flips/prove_*.out`, `depthcheck_*.out`, `nondet_merge.out`, `repro_merge.out`
- Production tables under comparison: `process/AGE_SY_v3/RefAlt.chrX.txt` (non-tiled `-t`) vs `process/AGE_SY_v3_tiled/RefAlt.chrX.txt` (tiled `-r`), both `-d 1000`

## The 7 flip sites

`16025382 16030427 16030958 16031180 16031851 16045058 16363526`

---

## Scripts, in run order — question → answer

1. **factorial_t_r_cap.chrX.sh** — `-t` vs `-r` over the *same* 500 kb region.
   **ANSWER: `-t` == `-r` identical.** The FLAG is not the cause. (Caveat: it gave
   both the same window, so it could not see the tile-vs-whole-chr *extent* effect.)

2. **test_3rd_allele_cap.chrX.sh** — hypothesis: rare 3rd allele + the cap.
   Sets up the diff of the two production RefAlt tables.

3. **allele_counts_at_flips.chrX.sh** — per-sample counts at the 7 sites.
   **SUPERSEDED / unfaithful:** used `-R` *isolated single bases* (zero flanking
   context), so its verdicts do not match the pipeline. Do not trust its kept/dropped.

4. **prove_depth_flips.chrX.sh** — EXACT caller, 500 kb window
   `chrX:15900000-16400000`, d=1000 vs d=50000, classified each flip.
   **ANSWER: 5/7 flips are the QUAL>59 gate, 1 is the 3rd allele, 1 is no-call.**
   Per-site (500 kb window):

   | site | d1000 | d50000 | gate |
   |---|---|---|---|
   | 16025382 | 63.8933 KEPT | 45.6247 DROP | B (QUAL) |
   | 16030427 | no-call | 180.213 KEEP | C (no-call) |
   | 16030958 | nALT2 635.87 DROP | nALT1 320.263 KEEP | A (3rd allele) |
   | 16031180 | 50.7705 DROP | 99.064 KEEP | B |
   | 16031851 | 56.7382 DROP | 123.154 KEEP | B |
   | 16045058 | 64.0468 KEEP | 58.0733 DROP | B |
   | 16363526 | 49.6448 DROP | 82.6348 KEEP | B |

5. **depth_cap_check.chrX.sh** — ±10 kb windows
   (`chrX:16015382-16055058,chrX:16353526-16373526`), true (uncapped, samtools
   -d 100000) per-sample depth at each site.
   **ANSWER: cap binds at 4 sites, NOT at 3.**

   | site | max true depth | samples ≥1000 (capped) |
   |---|---|---|
   | 16025382 | 1029 | 1 |
   | 16030427 | 1400 | 3 |
   | 16030958 | 2095 | 8 |
   | 16031180 | 1106 | 2 |
   | 16031851 | 566 | **0** |
   | 16045058 | 271 | **0** |
   | 16363526 | 759 | **0** |

   QUAL in THIS (±10 kb) window — **differs from the 500 kb window above**, which
   proves QUAL depends on the calling window even mid-window:

   | site | d1000 | d50000 |
   |---|---|---|
   | 16025382 | 36.3391 | 66.206 |
   | 16030427 | 96.212 | 132.212 |
   | 16030958 | 390.393 (nALT1) | 376.591 (nALT2) |
   | 16031180 | 110.727 | 85.5538 |
   | 16031851 | 134.905 | 111.443 |
   | 16045058 | 64.0468 | 67.0423 |
   | 16363526 | 24.2829 | 36.8386 |

6. **reproducibility_check.chrX.sh** — same ±10 kb window, d=1000 ×2 and
   d=50000 ×2 (parallel array + merge). Expected values baked in from #5.
   **STATUS: built, not yet the definitive run** (superseded in intent by #7).

7c. **account_59.chrX.sh** — per-SNP accounting of all 59 disagreements, NO rerun.
   KEY: the whole-chr caller KEEPS its BCF (bam2bcf2REFALT.sh writes calls.chrX.bcf,
   never deletes) — it's PRE the QUAL>59 filter, so it holds every called variant
   incl. those with QUAL<=59 and 3rd alleles (nALT>=2). The tiled caller DELETES its
   per-tile BCF (rm -f $tmpbcf). So the accounting is ONE-WAY: for the ~29 NEW-only
   SNPs (tiled kept, whole-chr dropped) the whole-chr BCF says exactly why whole-chr
   dropped it — 3rd allele (nALT>=2 -> -m2-M2), or biallelic QUAL<=59, or no-call.
   Near-cap sample (ref+alt >= ~900, ceiling ~940) = subsampling was in play; no
   capped sample = residual (identical reads yet flips = window effect on QUAL).
   OLD-only rows are inferential only (tiled BCF gone). **STATUS: RUN — see below.**
   The user accepted a one-way explanation (~30 SNPs) as sufficient.

   RESULT (account_59_54343282.out), 59 disagreements accounted:
   - tally: 27 RESIDUAL, 17 EXPLAINED-cap (OLD-only, inferential), 12 EXPLAINED-QUAL
     (NEW-only), 2 NO-CALL, 1 EXPLAINED-3rd.
   - EVERY SNP is mid-tile: min dL=1,018,399, min dR=85,871 — none within the 10kb
     pad. So NOT an edge/padding effect. (Confirms the user's original observation.)
   - Rigorous half = 29 NEW-only (whole-chr BCF is the dropping side):
     ~12 near-cap (wmaxDP 905-981, always the R8/Con_R8 deepest libraries) ->
       subsampling shifted whole-chr QUAL <=59 (e.g. 23207874 QUAL 3.5).
     ~15 UNCAPPED (wmaxDP down to 260, reads_same=YES) -> identical read counts yet
       whole-chr QUAL<=59 while tiled kept it -> window moves QUAL on identical data.
       Not hairline: 14279293 whole-chr QUAL=10.5, tiled kept (>59) — ~48pt swing.
     + 1 uncapped 3rd-allele (8868954, window changed the allele call), 2 no-call.
   - So ~half the disagreements = -d 1000 cap on the deepest libraries; ~half =
     bcftools QUAL depending on the calling window with identical reads (mechanism
     still unexplained). 59/250,830 = 0.02%.
   - Caveats: OLD-only EXPLAINED-cap are inferential (tiled BCF deleted); reads_same
     compares max+total depth, not the full 80-sample vector; wmaxDP 900-940 is
     ambiguous capped-vs-not (needs true uncapped depth to resolve).
   NEXT (optional): (1) airtight residual = full per-sample AD compare + samtools
   true depth at the 59, one SLURM job; (2) FIX = count-based SNP filter (counts are
   identical across callers; only pooled QUAL is fragile).

8. **forensic_residual.chrX.sh** — WHY does QUAL differ on identical reads? The
   residual "window changes QUAL" is a restatement, not a cause. QUAL comes from
   genotype likelihoods -> base qualities -> BAQ (base alignment quality, ON by
   default, no -B). Hypothesis: BAQ is computed window-dependently. Calls 16045058
   at window half-widths {50k,250k,1M,2.5M} x {BAQ on, -B off}, dumps QUAL + full
   per-sample AD each time. Reads: ADeq=YES across all => data identical; QUAL
   varies down BAQ-on col => window moves QUAL on same data; BAQ-off col FLAT =>
   BAQ is the cause. **STATUS: RUN (task 7=2.5M-BAQon TIMEOUT@8h, rest done; 5Mb
   arm dropped anyway). ANSWER: BAQ IS THE CAUSE.**

   RESULT (16045058, per-sample AD md5 + QUAL):
     condition     QUAL      nALT  AD_md5
     BAQon_50k     53.1493   1     8104e8c3
     BAQoff_50k    291.812   1     525d29f8
     BAQon_250k    64.0468   1     8104e8c3
     BAQoff_250k   294.81    1     525d29f8
     BAQon_1M      64.0468   1     8104e8c3
     BAQoff_1M     294.81    1     525d29f8
     BAQoff_2.5M   294.81    1     525d29f8
   - AD md5 = ONE value across ALL window sizes within a BAQ setting => reads are
     byte-identical regardless of window. Window size does NOT change the data.
   - BAQ OFF: QUAL ~292-295, flat across windows.
   - BAQ ON : QUAL 53->64, swings ACROSS the 59 cutoff (53.15@50k, 64.05@250k/1M).
   => the residual flip is a pure BAQ artifact: BAQ (local realignment, on by
      default, no -B) recomputes base qualities window-dependently; those feed
      PL/QUAL. Identical reads, different BAQ-adjusted base quals, QUAL crosses 59.
   - BAQ suppresses QUAL ~5x here (294 -> ~58): this site sits in context BAQ
      distrusts (low-complexity/near-indel), and QUAL>59 cuts through the wiggle.
   - Gap (small): BAQ-on QUAL stable by 250k (250k==1M==64.0468) but 50k=53.15 and
      whole-chr -t=58.07; a small method/small-window residual remains, all within
      the BAQ-on regime. Big picture settled: NO BAQ, NO FLIP.

   ============================ BOTTOM LINE ============================
   The tiled-vs-whole-chr SNP disagreements (59/250,830 = 0.02%, ALL mid-tile)
   have TWO causes, both acting on the pooled QUAL>59 filter (never on the counts,
   which are identical everywhere):
     (1) ~half: -d 1000 cap on the deepest libraries (R8/Con_R8) -> different read
         draw per window -> QUAL shifts.
     (2) ~half: BAQ recomputes base qualities window-dependently on IDENTICAL reads
         -> QUAL shifts (this experiment).
   FIX: select SNPs on the depth/window-stable allele counts, not pooled QUAL>59.
   (Optionally add -B to disable BAQ, but that changes calls globally.)

   Result so far (maxdepth_at_diff_snps): 250,771 SNPs agree, 59 disagree. The
   disagreements sit at high depth: ~31/59 have a sample >=900 (near the ~940
   capped ceiling) => subsampling-plausible; ~28/59 have max per-sample depth <900
   => NO capped sample => residual window/QUAL effect (e.g. 16045058 AB8, true depth
   271). So it is NOT all sampling — roughly half cap, half residual. (My first
   read said "7 residual" using a wrong <500 threshold; the capped value in these
   -d1000 tables is ~940, so <900 => uncapped, giving ~28.)

7b. **maxdepth_at_diff_snps.chrX.sh** — the RIGHT count-based analysis, NO rerun.
   For the actual tiled-vs-whole disagreeing SNPs (from the two existing RefAlt
   tables), computes MAX per-sample depth (ref+alt) and asks: is a sample sitting
   near the -d 1000 ceiling? Disagreeing SNPs with max < 500 have no plausibly
   capped sample → both callers saw identical reads → the disagreement is NOT
   subsampling (residual). **STATUS: built, not yet run. Run on login node.**
   This is the honest test of "can we explain every difference by sampling?"
   NOTE: characterize_diff_snps.sh tracked the WRONG summary (summed reads +
   alt-freq); the cap is per-sample so MAX per-sample depth is what matters.

7. **nondeterminism_check.chrX.sh** — the EXACT same command run TWICE (500 kb
   window, same `-d`) and diff the FULL filtered SNP tables (d1000 A/B, d50000 A/B).
   **ANSWER: DETERMINISTIC.** `nondet_merge.out`: d=1000 → 6723 vs 6723 SNPs, 0 diff;
   d=50000 → 6726 vs 6726 SNPs, 0 diff. Identical command → byte-identical result.
   So the tiled/non-tiled disagreement is NOT run-to-run noise — it is a real,
   reproducible **deterministic window-extent effect** (candidate b).

---

## Conclusions so far (established, not speculation)

- Dropped SNPs are dropped by the **`QUAL>59` gate** (5/7) or the **`-m2 -M2`
  biallelic gate** when a rare 3rd allele appears/vanishes (1/7). Immediate cause = the filter.
- **Per-sample allele counts (AD) are stable** across depth; **pooled QUAL is not.**
- **QUAL moves with `-d` AND with the calling window**, including:
  - at sites where **no sample hits the cap** (16045058: max depth 271, uncapped at
    both `-d`, yet QUAL 64.05 vs 67.04). So `-d` perturbs QUAL through a path that is
    **not** read-subsampling. Mechanism unknown — do NOT guess.
  - **mid-window** (16031851 sits ~130 kb inside the 500 kb window; QUAL 56.7 there
    vs 134.9 in the ±10 kb window).
- **The per-sample AD counts are stable across the WINDOW too, not just across `-d`.**
  Same sample, same site, 500 kb window vs ±10 kb window: AD identical
  (16030958 AgeSY10_R8_F d1000 = 944:711,228 in both; 16030427 AgeSY10_R8_F d50000
  = 1170:1166,0 in both). Only pooled QUAL moves with the window; the counts do not.
- **REJECTED: boundary/padding effect.** The tiled/non-tiled callers differ at MANY
  SNPs mid-tile, nowhere near edges — padding only affects near-edge positions.
- **REJECTED: nondeterminism.** Identical command twice = byte-identical (see #7).

## Open question — RESOLVED to "deterministic window effect"

Was: (a) nondeterministic caller vs (b) deterministic global window-extent effect.
`nondeterminism_check` (#7) answered **(b)**: identical command → identical result,
so the tiled/non-tiled disagreement is a real, reproducible consequence of the
window each caller uses (5 Mb tile vs whole-chr), NOT noise. The internal mechanism
(why window extent changes a mid-window, uncapped QUAL) is still not explained — but
it does not need to be, because the fix below does not depend on it.

## Practical fix (holds regardless of the open question)

The RefAlt allele counts are stable across BOTH `-d` and window; only the pooled
`QUAL>59` (and the biallelic `-m2 -M2`) filter is depth/window-fragile. The durable
fix: select SNPs on a **count-based, depth/window-robust** criterion instead of
pooled QUAL. NOT YET BUILT.

**Validation is cheap and needs NO caller re-run:** apply the count-based rule to
the two EXISTING production tables `AGE_SY_v3/RefAlt.chrX.txt` (non-tiled) and
`AGE_SY_v3_tiled/RefAlt.chrX.txt` (tiled). If the rule makes their SNP sets agree,
the fix works and the tiled/non-tiled discrepancy is closed.
