#!/bin/bash
#SBATCH --job-name=filter_grid
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:20:00
#SBATCH -o logs/AGE_SY/filter_grid.out
# How many SNPs survive for each (indel-distance x maxaf), holding >=10 reads and
# differs-between-founders fixed. Pure count off catalog.annot.tsv.gz, no recount.
# B5 ignored on chr2L.  ->  logs/AGE_SY/filter_grid.out
set -uo pipefail
ANNOT=process/AGE_SY_v5/Catalog/catalog.annot.tsv.gz
[[ -s "$ANNOT" ]] || { echo "missing $ANNOT" >&2; exit 1; }

zcat "$ANNOT" | awk -F'\t' '
BEGIN{
  ng=split("0 5 10 15 20 25 30 40 50",G," ");
  nm=split("0.03 0.05 0.07",M," ");
}
NR==1{ for(i=6;i<=NF;i++){n=$i; sub(/^AD_/,"",n); if(n=="B5") b5=i} next }
{
  chr=$1; d=$5+0; cand++;
  depthok=1;
  for(k=1;k<=nm;k++){ nf[k]=1; hr[k]=0; ha[k]=0 }
  for(i=6;i<=NF;i++){
    if(chr=="chr2L" && i==b5) continue;             # B5 exempt on chr2L
    split($i,a,","); dep=a[1]+a[2];
    if(dep<10) depthok=0;
    af=(dep>0? a[2]/dep : -1);
    for(k=1;k<=nm;k++){ m=M[k];
      if(af<=m) hr[k]=1; else if(af>=1-m) ha[k]=1; else nf[k]=0 }
  }
  if(!depthok) next;
  for(k=1;k<=nm;k++){
    if(!(nf[k] && hr[k] && ha[k])) continue;         # near-fixed + differs
    for(g=1;g<=ng;g++) if(d>=G[g]) cnt[g,k]++;        # indel distance
  }
}
END{
  printf "candidate SNPs (biallelic in founders): %d\n\n", cand;
  printf "SNPs kept  (rows = >= bp from indel, cols = maxaf; >=10 reads + differs held)\n";
  printf "%-10s", "indel>=";
  for(k=1;k<=nm;k++) printf "%12s", "af"M[k];
  print "";
  for(g=1;g<=ng;g++){ printf "%-10s", G[g];
    for(k=1;k<=nm;k++) printf "%12d", cnt[g,k];
    print "" }
  print "\n(current default = indel>=25, af0.03)";
}'
