#!/bin/bash
#SBATCH --job-name=flag_ladder
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=02:00:00
#SBATCH -o logs/AGE_SY/flag_ladder_%j.out
###############################################################################
# flag_ladder.chrX.sh — attribute the v3-vs-v4 ALT loss to a specific mpileup flag.
#
# ONE sample, chrX catalog only. Re-counts under a ladder of mpileup settings, all
# with the SAME caller (call -m -C alleles -T catalog), so only the flag changes.
# Reports summed REF/ALT over the chrX catalog sites for each setting.
#
#   v4       : -B -q20 -Q20      (current catalog caller)
#   -q off   : -B -q0  -Q20      (isolates mapping-quality filter)
#   -Q off   : -B -q20 -Q13      (isolates base-quality filter)
#   pivot    : -B -q0  -Q13      (both filters off, BAQ still off)
#   v3-like  :    -q0  -Q13      (pivot + BAQ ON = v3 filter levels; pivot-vs-this = BAQ)
#
#   sbatch scripts_oneoffs/AGE_SY/catalog_v4/flag_ladder.chrX.sh [sample_bam]
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true
REF=pipeline/ref/dm6.fa
CAT=process/AGE_SY_v4/Catalog/catalog.tsv.gz
BAM=${1:-data/bam/AGE_SY/Con_R5_F.bam}     # top-offender female pool (high chrX cov)
for f in "$REF" "$CAT" "$BAM"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done
echo "sample: $BAM   region: chrX   catalog: $CAT"

count() {   # $1=label  $2..=extra mpileup flags
  local lab="$1"; shift
  bcftools mpileup "$@" -r chrX -T "$CAT" -a FORMAT/AD -f "$REF" "$BAM" 2>/dev/null \
    | bcftools call -m -C alleles -T "$CAT" 2>/dev/null \
    | bcftools query -f '[%AD]\n' 2>/dev/null \
    | awk -F, 'BEGIN{r=0;a=0} {r+=$1; a+=$2} END{printf "  %-9s  REF=%-10d ALT=%-9d  ALT%%=%.3f\n", L, r, a, 100*a/(r+a)}' L="$lab"
}
echo "=== summed REF/ALT over chrX catalog sites (same caller, only flag varies) ==="
count "v4"       -B -q 20 -Q 20
count "-q off"   -B -q 0  -Q 20
count "-Q off"   -B -q 20 -Q 13
count "pivot"    -B -q 0  -Q 13
count "v3-like"     -q 0  -Q 13
echo
echo "read: v4->pivot moving = the -q/-Q filters; -q off vs -Q off = which filter;"
echo "      pivot vs v3-like = BAQ.  Compare ALT% to the real v3/v4 tables for this sample."
