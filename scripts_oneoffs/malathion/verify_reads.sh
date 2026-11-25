#!/bin/bash
#SBATCH --job-name=mal_verify
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=2:00:00
#SBATCH -o logs/malathion/verify_%j.out
#SBATCH -e logs/malathion/verify_%j.out

# verify_reads.sh — after the download on the transfer host, name the files the
# way fq2bam.sh expects and check every one against ENA's md5.
#
# The download is one command typed on access-hpc3.rcic.uci.edu (restricted
# shell, absolute paths, no cd):
#   wget -c -i /dfs7/adl/tdlong/fly_pool/XQTL2-dev/helpfiles/malathion/sra_urls.txt -P /dfs7/adl/tdlong/fly_pool/XQTL2-dev/data/raw/malathion_2022
# which leaves ENA's names, <run>_1.fastq.gz and <run>_2.fastq.gz. fq2bam.sh finds
# a sample by the barcode pair "<colA>-<colB>" and "-R1." in the filename, and
# the barcodes file gives <run> and SRA, so each file is renamed to
# <run>-SRA-R1.fastq.gz / -R2. Then md5 against helpfiles/malathion/sra_runs.tsv;
# a match gets a <file>.md5ok marker, a mismatch is deleted so wget -c cannot
# resume a corrupt file. rerun_2022_from_reads.sh needs all 16 markers.
#
# Submit from the repo root on the login node once the wget has finished:
#   sbatch scripts_oneoffs/malathion/verify_reads.sh
# Reads 29 GB, so it is a job, not a login-node command. Rerun as often as needed.

TSV=helpfiles/malathion/sra_runs.tsv
RAW=data/raw/malathion_2022
[ -f "$TSV" ] || { echo "run from the XQTL2-dev root ($TSV not found)" >&2; exit 1; }
[ -d "$RAW" ]  || { echo "no $RAW yet" >&2; exit 1; }

n_ok=0; n_bad=0; n_missing=0
while IFS=$'\t' read -r run sample alias url1 url2 md5_1 md5_2 bytes1 bytes2; do
  [ "$run" = "run" ] && continue
  for i in 1 2; do
    md5="md5_$i"; md5=${!md5}; bytes="bytes$i"; bytes=${!bytes}
    ena="$RAW/${run}_$i.fastq.gz"
    dest="$RAW/${run}-SRA-R$i.fastq.gz"
    [ -e "$ena" ] && [ ! -e "$dest" ] && mv "$ena" "$dest"
    if [ ! -s "$dest" ]; then
      echo "MISSING  $dest  ($sample)"; n_missing=$((n_missing+1)); continue
    fi
    if [ -e "$dest.md5ok" ]; then
      echo "ok       $dest  (verified earlier)"; n_ok=$((n_ok+1)); continue
    fi
    size=$(stat -c %s "$dest")
    if [ "$size" -ne "$bytes" ]; then
      echo "SHORT    $dest  $size of $bytes bytes -- incomplete; rerun the wget (it resumes)"; n_bad=$((n_bad+1)); continue
    fi
    got=$(md5sum "$dest" | cut -d' ' -f1)
    if [ "$got" = "$md5" ]; then
      touch "$dest.md5ok"; echo "ok       $dest  ($sample, $(du -h "$dest" | cut -f1))"; n_ok=$((n_ok+1))
    else
      echo "CORRUPT  $dest  md5 $got, expected $md5 -- deleted; rerun the wget"; rm -f "$dest"; n_bad=$((n_bad+1))
    fi
  done
done < "$TSV"

echo
echo "verified $n_ok / 16   incomplete or corrupt $n_bad   missing $n_missing"
if [ "$n_ok" -eq 16 ]; then
  echo "all 16 verified. next, from the repo root:  bash scripts_oneoffs/malathion/rerun_2022_from_reads.sh"
  exit 0
fi
echo "not complete. rerun the wget on access-hpc3 (it resumes), then: sbatch scripts_oneoffs/malathion/verify_reads.sh"
exit 1
