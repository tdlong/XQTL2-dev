AGE_SY Processing Notes
=======================

Project: AGE_SY (aging experiment with SY10 and SY20 treatments)
Groups: Con (control), AgeSY10, AgeSY20
Replicates: R1-R12 (both sexes: F, M)
Final BAMs land in: data/bam/AGE_SY/

BAM naming convention: {Group}_{Replicate}_{Sex}.bam
  e.g. Con_R8_F.bam, AgeSY10_R9_M.bam, AgeSY20_R11_F.bam


--- R1-R7: already processed ---

BAMs for R1-R7 (all groups, both sexes) already exist in data/bam/AGE_SY/.
No reprocessing needed.


--- R8-R11: two sequencing runs, merged ---

Samples R8-R11 were sequenced twice across two runs:
  March 2026:  helpfiles/AGE_SY/readname.mapping.AGE_SY.March_26.txt
               raw FASTQs: /dfs7/adl/tdlong/fly_pool/XQTL2-dev/data/raw/AGE_SY/March_26
  April 2026:  helpfiles/AGE_SY/readname.mapping.AGE_SY.April5_26.txt
               raw FASTQs: /dfs7/adl/tdlong/fly_pool/XQTL2-dev/data/raw/AGE_SY/April5_26

Step 1: align each run to a temporary BAM directory

  mkdir -p data/bam/AGE_SY_tempA
  JOBID_MARCH=$(sbatch --array=1-24 --parsable pipeline/scripts/fq2bam.sh \
      helpfiles/AGE_SY/readname.mapping.AGE_SY.March_26.txt \
      /dfs7/adl/tdlong/fly_pool/XQTL2-dev/data/raw/AGE_SY/March_26 \
      data/bam/AGE_SY_tempA)

  mkdir -p data/bam/AGE_SY_tempB
  JOBID_APRIL=$(sbatch --array=1-24 --parsable pipeline/scripts/fq2bam.sh \
      helpfiles/AGE_SY/readname.mapping.AGE_SY.April5_26.txt \
      /dfs7/adl/tdlong/fly_pool/XQTL2-dev/data/raw/AGE_SY/April5_26 \
      data/bam/AGE_SY_tempB)

Step 2: merge tempA + tempB into the main BAM directory

  sbatch --dependency=afterok:${JOBID_MARCH}:${JOBID_APRIL} \
      scripts_oneoffs/AGE_SY/round2_v2_R8-R11/mergebams.sh

  (Or submit manually once both arrays have completed.)
  mergebams.sh loops through tempA, finds the matching sample in tempB,
  and writes merged+indexed BAMs to data/bam/AGE_SY/.

Step 3: verify merged BAMs, then remove temp directories

  Once you are satisfied the merged BAMs in data/bam/AGE_SY/ are correct:

  rm -rf data/bam/AGE_SY_tempA data/bam/AGE_SY_tempB


--- Updating helpfiles after merge ---

Step 4: check BAM file sizes

  The pipeline excludes BAMs under 1 GB (low coverage causes downstream
  problems). Check sizes of the new R8-R11 BAMs on the cluster:

  ls -l data/bam/AGE_SY/{AgeSY10,AgeSY20,Con}_R{8,9,10,11}_{F,M}.bam

  Any BAM under 1 GB should be excluded from the analysis.

Step 5: update the four design files

  Design files are generated from helpfiles/AGE_SY/summary_info_v1.xlsx.
  The script scripts_oneoffs/AGE_SY/common/make_AGE_SY_design_files.py reads the
  Control_Flies and Aged_Flies tabs to get Num and Proportion values.

  Run from XQTL2-dev root (requires openpyxl):
    python3 scripts_oneoffs/AGE_SY/common/make_AGE_SY_design_files.py

  This writes filesize=0 as a placeholder. After getting BAM sizes from
  Step 4, replace the 0s with actual sizes. Also remove any rows for BAMs
  that failed the 1 GB threshold.

Step 6: update AGE_SY.bams

  Add R8-R11 BAM paths (only those passing the 1 GB threshold) to
  helpfiles/AGE_SY/AGE_SY.bams. The file already contains R1-R7 sample
  BAMs and the B-population founder BAMs at the bottom — insert the new
  sample paths before the founder lines.

Step 7: update AGE_SY_haplotype_parameters.R

  Add passing R8-R11 sample names to names_in_bam. The pipeline README
  gives a one-liner to regenerate names_in_bam from the BAM directory
  (filters by -size +1G automatically):

    echo -n "names_in_bam <- c(" && \
    find data/bam/AGE_SY -name "*.bam" -size +1G -print0 | \
    xargs -0 -n1 basename | sed 's/.bam//' | sort | \
    sed 's/.*/"&"/' | tr '\n' ',' | sed 's/,$//' && echo ")"

  Paste the output into AGE_SY_haplotype_parameters.R replacing the
  existing names_in_bam line.

  Commit and push all updated helpfiles before submitting the pipeline.


--- Running the pipeline (from REFALT onward) ---

  Prior results (R1-R7 only) are in process/AGE_SY/.
  This run (R1-R11) uses process/AGE_SY_v2/ to preserve old results.

  mkdir -p process/AGE_SY_v2

  JID_REFALT=$(bash pipeline/scripts/run_refalt.sh \
      --bamlist helpfiles/AGE_SY/AGE_SY.bams \
      --dir     process/AGE_SY_v2)

  JID_HAPS=$(bash pipeline/scripts/run_haps.sh \
      --after   $JID_REFALT \
      --parfile helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R \
      --dir     process/AGE_SY_v2)

  bash pipeline/scripts/run_scan.sh \
      --after $JID_HAPS --dir process/AGE_SY_v2 \
      --scan AGE_SY10_F --design helpfiles/AGE_SY/AGE_SY10_F.test.txt

  bash pipeline/scripts/run_scan.sh \
      --after $JID_HAPS --dir process/AGE_SY_v2 \
      --scan AGE_SY10_M --design helpfiles/AGE_SY/AGE_SY10_M.test.txt

  bash pipeline/scripts/run_scan.sh \
      --after $JID_HAPS --dir process/AGE_SY_v2 \
      --scan AGE_SY20_F --design helpfiles/AGE_SY/AGE_SY20_F.test.txt

  bash pipeline/scripts/run_scan.sh \
      --after $JID_HAPS --dir process/AGE_SY_v2 \
      --scan AGE_SY20_M --design helpfiles/AGE_SY/AGE_SY20_M.test.txt


===========================================================================
--- R12: round 3, single run (July_26) ---
===========================================================================

Scripts:  scripts_oneoffs/AGE_SY/round3_v3_R12/
Results:  process/AGE_SY_v3/   (round 1 = process/AGE_SY [R1-R7],
                                round 2 = process/AGE_SY_v2 [R1-R11])

R12 was a SINGLE sequencing run (xR096, July_26), so unlike R8-R11 there is
NO two-run merge — fq2bam writes straight into data/bam/AGE_SY/.

Barcode / readname mapping (one file per run, in helpfiles/AGE_SY/):
  readname.mapping.AGE_SY.July_26.txt   (<i7>  <i5>  <sample>; 6 rows)

Step 1: align (single run, no merge)

  bash scripts_oneoffs/AGE_SY/round3_v3_R12/01_align.sh
  # = sbatch --array=1-6 pipeline/scripts/fq2bam.sh \
  #       helpfiles/AGE_SY/readname.mapping.AGE_SY.July_26.txt \
  #       /dfs7/adl/tdlong/fly_pool/XQTL2-dev/data/raw/AGE_SY/July_26 \
  #       data/bam/AGE_SY

Step 2: update helpfiles for R12 (same as Steps 4-7 above, now for R12)
  - check R12 BAM sizes (>1 GB) in data/bam/AGE_SY/
  - add R12 fly counts to summary_info_v1.xlsx, rerun make_AGE_SY_design_files.py
    (REPS already bumped to R1-R12)
  - add passing R12 BAMs to AGE_SY.bams; regenerate names_in_bam
  - commit + push, then git pull on the cluster

Step 3: recall SNPs -> haplotypes -> scan (into process/AGE_SY_v3)

  bash scripts_oneoffs/AGE_SY/round3_v3_R12/02_refalt_haps_scan.sh
  # REFALT -> haps -> the four scans (AGE_SY10/20 x F/M), same as the v2 block above
