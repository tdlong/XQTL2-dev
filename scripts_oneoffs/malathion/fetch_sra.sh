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

# fetch_sra.sh — the eight malathion pools of Macdonald & Long 2022 (BioProject
# PRJNA857080) from the ENA mirror, with wget, each file checked against ENA's
# md5. No sra-tools needed.
#
# Runs as a batch array, one run per task, both mates. A task exits non-zero
# unless both files match ENA's md5, so a chain that waits on it with afterok
# never starts on a partial download. Submitted by rerun_2022_from_reads.sh;
# by hand from the repo root:
#   sbatch scripts_oneoffs/malathion/fetch_sra.sh
# Without SLURM (plain `bash`, e.g. on a transfer host) it fetches all eight in
# parallel and logs to logs/malathion/fetch_dtn.out.
#
# Files land in data/raw/malathion_2022 as <run>-SRA-R1.fastq.gz and -R2, which is
# the naming fq2bam.sh looks for: it greps the barcode pair "<colA>-<colB>" and
# "-R1." out of the filename, and helpfiles/malathion/malathion_2022.barcodes.txt
# gives <run> and SRA as the two "barcodes". A verified file gets a <file>.md5ok
# marker; rerun_2022_from_reads.sh skips the fetch when all sixteen are there.
# Idempotent: a file with a marker, or with the right md5, is not fetched again.
# Never run two copies at once: two wgets into one file corrupt it.

TSV=helpfiles/malathion/sra_runs.tsv
OUT=data/raw/malathion_2022
[ -f "$TSV" ] || { echo "run from the XQTL2-dev root ($TSV not found)" >&2; exit 1; }
mkdir -p "$OUT" logs/malathion

fetch_one () {  # fetch_one <row number 1-8>; returns 1 on any failure
  local line run sample alias url1 url2 md5_1 md5_2 bytes1 bytes2 i url md5 dest got fail=0
  line=$(sed -n "$(($1 + 1))p" "$TSV")              # line 1 is the header
  [ -n "$line" ] || { echo "no row $1 in $TSV" >&2; return 1; }
  IFS=$'\t' read -r run sample alias url1 url2 md5_1 md5_2 bytes1 bytes2 <<< "$line"
  echo "$(date '+%F %T')  $run = $sample ($alias)"
  for i in 1 2; do
    url="url$i"; md5="md5_$i"; url=${!url}; md5=${!md5}
    dest="$OUT/${run}-SRA-R$i.fastq.gz"
    if [ -e "$dest.md5ok" ] && [ -s "$dest" ]; then
      echo "  have $dest (verified earlier)"; continue
    fi
    if [ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5" ]; then
      echo "  have $dest (md5 ok)"; touch "$dest.md5ok"; continue
    fi
    echo "  fetching https://$url -> $dest"
    wget -c -q -O "$dest" "https://$url" || { echo "  wget failed for $url" >&2; fail=1; continue; }
    got=$(md5sum "$dest" | cut -d' ' -f1)
    if [ "$got" = "$md5" ]; then
      touch "$dest.md5ok"
      echo "  ok   $dest  $(du -h "$dest" | cut -f1)"
    else
      echo "  MD5 MISMATCH $dest: got $got want $md5 -- removing" >&2
      rm -f "$dest" "$dest.md5ok"; fail=1
    fi
  done
  return $fail
}

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  # one run per array task
  fetch_one "$SLURM_ARRAY_TASK_ID"; rc=$?
  echo "$(date '+%F %T')  done, exit $rc"
  exit $rc
fi

# plain shell on the data-transfer host: all eight in parallel, one log
exec >> logs/malathion/fetch_dtn.out 2>&1
echo "=== $(date '+%F %T')  $(hostname)  fetching 8 runs of PRJNA857080 ==="
pids=(); rows=()
for t in 1 2 3 4 5 6 7 8; do fetch_one "$t" & pids+=($!); rows+=($t); done
fail=0
for k in "${!pids[@]}"; do wait "${pids[$k]}" || { echo "row ${rows[$k]} FAILED"; fail=1; }; done
n_ok=$(ls "$OUT"/*.md5ok 2>/dev/null | wc -l)
echo "=== $(date '+%F %T')  done: $n_ok / 16 files verified, exit $fail ==="
exit $fail
