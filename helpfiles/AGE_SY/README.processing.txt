AGE_SY Processing Notes
=======================

Project: AGE_SY (aging experiment with SY10 and SY20 treatments)
Groups: Con (control), AgeSY10, AgeSY20
Replicates: R1-R11 (both sexes: F, M)
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
      scripts_oneoffs/mergebams.sh

  (Or submit manually once both arrays have completed.)
  mergebams.sh loops through tempA, finds the matching sample in tempB,
  and writes merged+indexed BAMs to data/bam/AGE_SY/.

Step 3: verify merged BAMs, then remove temp directories

  Once you are satisfied the merged BAMs in data/bam/AGE_SY/ are correct:

  rm -rf data/bam/AGE_SY_tempA data/bam/AGE_SY_tempB
