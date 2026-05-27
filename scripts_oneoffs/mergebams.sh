#!/bin/bash
#SBATCH --job-name=merge_bams
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

module load samtools/1.10

TEMP_A="data/bam/AGE_SY_tempA"
TEMP_B="data/bam/AGE_SY_tempB"
OUTPUT_DIR="data/bam/AGE_SY"

SAMPLE=$(ls ${TEMP_A}/*.bam | sed -n "${SLURM_ARRAY_TASK_ID}p" | xargs basename | sed 's/.bam//')

BAM_A="${TEMP_A}/${SAMPLE}.bam"
BAM_B="${TEMP_B}/${SAMPLE}.bam"
OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE}.bam"

echo "Task ${SLURM_ARRAY_TASK_ID}: merging ${SAMPLE}"

if [[ ! -f ${BAM_B} ]]; then
    echo "ERROR: ${BAM_B} not found" >&2
    exit 1
fi

samtools merge -f ${OUTPUT_BAM} ${BAM_A} ${BAM_B}
samtools index ${OUTPUT_BAM}

echo "Done: ${OUTPUT_BAM}"
