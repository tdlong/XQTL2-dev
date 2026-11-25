#!/bin/bash
#SBATCH --job-name=sra_fetch
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
#SBATCH --array=1-8
#SBATCH -o logs/malathion/fetch_%A_%a.out
#SBATCH -e logs/malathion/fetch_%A_%a.out

# fetch_sra.sh — one array task per run of BioProject PRJNA857080 (the eight
# malathion pools of Macdonald & Long 2022), pulled from the ENA mirror with
# wget and checked against ENA's md5. No sra-tools needed.
#
# Files land in data/raw/malathion_2022 as <run>-SRA-R1.fastq.gz and -R2, which is
# the naming fq2bam.sh looks for: it greps the barcode pair "<colA>-<colB>" and
# "-R1." out of the filename, and helpfiles/malathion/malathion_2022.barcodes.txt
# gives <run> and SRA as the two "barcodes". Idempotent: a file already present
# with the right md5 is not fetched again.
#
# Submitted by rerun_2022_from_reads.sh; by hand from the repo root:
#   sbatch scripts_oneoffs/malathion/fetch_sra.sh

TSV=helpfiles/malathion/sra_runs.tsv
OUT=data/raw/malathion_2022
mkdir -p "$OUT"

line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TSV")      # line 1 is the header
[ -n "$line" ] || { echo "no row $SLURM_ARRAY_TASK_ID in $TSV" >&2; exit 1; }
IFS=$'\t' read -r run sample alias url1 url2 md5_1 md5_2 bytes1 bytes2 <<< "$line"
echo "$(date)  task $SLURM_ARRAY_TASK_ID: $run = $sample ($alias)"

fail=0
for i in 1 2; do
  url="url$i"; md5="md5_$i"; url=${!url}; md5=${!md5}
  dest="$OUT/${run}-SRA-R$i.fastq.gz"
  if [ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5" ]; then
    echo "  have $dest (md5 ok)"; continue
  fi
  echo "  fetching https://$url -> $dest"
  wget -c -q -O "$dest" "https://$url" || { echo "  wget failed for $url" >&2; fail=1; continue; }
  got=$(md5sum "$dest" | cut -d' ' -f1)
  if [ "$got" = "$md5" ]; then
    echo "  ok   $dest  $(du -h "$dest" | cut -f1)"
  else
    echo "  MD5 MISMATCH $dest: got $got want $md5 -- removing" >&2
    rm -f "$dest"; fail=1
  fi
done
echo "$(date)  done, exit $fail"
exit "$fail"
