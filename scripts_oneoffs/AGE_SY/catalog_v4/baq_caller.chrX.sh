#!/bin/bash
#SBATCH --job-name=baq_caller
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:40:00
#SBATCH -o logs/AGE_SY/baq_caller_%a.out
###############################################################################
# baq_caller.chrX.sh — which setting causes v4's ALT loss: BAQ or the -C alleles
# constraint? (-q/-Q already shown irrelevant.) 2x2, ONE sample, chrX, correct
# --max-depth 2000, run as a PARALLEL 4-task array + dependent merge.
#
#   1 catalog : -B  call -m -C alleles   (BAQ off, constrained = real v4 recipe)
#   2 +BAQ    :     call -m -C alleles   (BAQ on,  constrained)   catalog->this = BAQ
#   3 +free   : -B  call -m               (BAQ off, unconstrained) catalog->this = -C constraint
#   4 v3ish   :     call -m               (BAQ on,  unconstrained ~ v3)
#
# This is a plain SLURM array script. Submit it yourself so the sbatch is visible:
#   JID=$(sbatch --parsable --array=1-4 scripts_oneoffs/AGE_SY/catalog_v4/baq_caller.chrX.sh)
#   sbatch --dependency=afterok:$JID scripts_oneoffs/AGE_SY/catalog_v4/baq_caller.chrX.sh --merge
# The array (1-4) runs the 4 settings in parallel; the --merge task combines them.
# BAM=... via --export to override the sample.
###############################################################################
set -uo pipefail
REF=pipeline/ref/dm6.fa
CAT=process/AGE_SY_v4/Catalog/catalog.tsv.gz
BAM=${BAM:-data/bam/AGE_SY/Con_R5_F.bam}
SAMP=$(basename "$BAM" .bam)
WORK=process/AGE_SY_v4/baq_caller_work
OUT=process/AGE_SY_v4/baq_caller.${SAMP}.chrX.tsv.gz
S1="-B|-C alleles|catalog"; S2="|-C alleles|baqon"; S3="-B||free"; S4="||v3ish"
mkdir -p "$WORK"

# -------- merge task ( sbatch ... baq_caller.chrX.sh --merge ) --------
if [[ "${1:-}" == "--merge" ]]; then
  module load bcftools/1.21 2>/dev/null || true
  { echo -e "POS\tref.catalog\talt.catalog\tref.baqon\talt.baqon\tref.free\talt.free\tref.v3ish\talt.v3ish"
    paste "$WORK/pos" "$WORK/1" "$WORK/2" "$WORK/3" "$WORK/4"; } | bgzip > "$OUT"
  echo "rows: $(zcat "$OUT" | wc -l)  -> $OUT"
  zcat "$OUT" | awk -F'\t' 'NR>1{cr+=$2;ca+=$3;br+=$4;ba+=$5;fr+=$6;fa+=$7;vr+=$8;va+=$9}
    END{printf "ALT%%  catalog(-B,-C)=%.3f  +BAQ(-C)=%.3f  +free(-B,-m)=%.3f  v3ish(-m,BAQon)=%.3f\n",
        100*ca/(cr+ca),100*ba/(br+ba),100*fa/(fr+fa),100*va/(vr+va)}'
  echo "read: catalog->+BAQ = BAQ effect;  catalog->+free = the -C alleles constraint;  v3ish ~ real v3"
  exit 0
fi

# -------- array task: one setting --------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: this is a SLURM array script — do not 'bash' it. Submit with:" >&2
  echo "  JID=\$(sbatch --parsable --array=1-4 $0); echo \"array \$JID\"" >&2
  echo "  sbatch --dependency=afterok:\$JID -o logs/AGE_SY/baq_caller_merge.out $0 --merge" >&2
  exit 1
fi
module load bcftools/1.21 2>/dev/null || true
eval "s=\$S$SLURM_ARRAY_TASK_ID"; IFS='|' read -r BAQ CONS LAB <<< "$s"
[[ "$SLURM_ARRAY_TASK_ID" == "1" ]] && \
  bcftools mpileup -B -q0 -Q13 --max-depth 2000 -r chrX -T "$CAT" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
    | bcftools call -m -C alleles -T "$CAT" 2>/dev/null | bcftools query -f '%POS\n' > "$WORK/pos"
echo "task $SLURM_ARRAY_TASK_ID ($LAB): BAQ='${BAQ:-on}' call -m $CONS" >&2
bcftools mpileup $BAQ -q 0 -Q 13 --max-depth 2000 -r chrX -T "$CAT" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
  | bcftools call -m $CONS -T "$CAT" 2>/dev/null \
  | bcftools query -f '[%AD]\n' \
  | awk -F, '{n=NF; r=($1=="."?0:$1); t=0; for(i=1;i<=n;i++) t+=($i=="."?0:$i); print r"\t"(t-r)}' > "$WORK/$SLURM_ARRAY_TASK_ID"
echo "task $SLURM_ARRAY_TASK_ID done"
