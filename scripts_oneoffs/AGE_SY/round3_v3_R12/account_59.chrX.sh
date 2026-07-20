#!/usr/bin/env bash
# account_59.chrX.sh  —  one row per disagreeing SNP, sorted into "explainable"
# and "residual", from data already on disk. No mpileup.
#
# The two mechanisms that CAN flip a call between the tiled and whole-chr runs
# both require a sample sitting near the -d 1000 ceiling (subsampling only bites
# there):
#   (a) 3rd allele + capped sample -> subsampling flips biallelic<->triallelic
#       (the biallelic filter -m2 -M2 then keeps one run, drops the other).
#   (b) rare/intermediate allele + capped sample -> subsampling nudges the pooled
#       QUAL across 59.
# So NEAR-CAP is the pivot. In the FILTERED tables a 3rd allele is invisible
# (all biallelic), but the whole-chr BCF (calls.chrX.bcf, KEPT, pre-filter) shows
# it directly as nALT>=2 -- for the SNPs whole-chr dropped. For the other
# direction the tiled BCF is gone, so there the smoking gun is a near-cap sample.
#
# Per SNP: caller, tile-edge distances, whole-chr QUAL/nALT, max per-sample depth
# (near ~940 => a sample was capped), reads-identical check, and a verdict:
#   EXPLAINED-3rd   whole-chr emitted a 3rd allele (nALT>=2)
#   EXPLAINED-QUAL  near-cap sample, biallelic, QUAL near 59 -> subsampling/QUAL
#   RESIDUAL        NO capped sample -> identical reads yet flips -> window/QUAL,
#                   mechanism unexplained (this is the honest leftover)
#
# Reads only (+ index the bcf once); login node:
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/account_59.chrX.sh
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true

A=process/AGE_SY_v3/RefAlt.chrX.txt          # OLD = whole-chr (-t)
B=process/AGE_SY_v3_tiled/RefAlt.chrX.txt     # NEW = tiled (-r)
BCF=process/AGE_SY_v3/calls.chrX.bcf          # whole-chr BCF (kept, PRE QUAL filter)
TILE=5000000
NEARCAP=900                                    # ref+alt >= this => a sample near the -d 1000 ceiling (~940)
for f in "$A" "$B" "$BCF"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done
[[ -e "$BCF.csi" || -e "$BCF.tbi" ]] || bcftools index "$BCF"

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# 1. disagreeing positions, tagged
awk 'NR>1{print $2}' "$A" | sort -u > "$TMP/a.pos"
awk 'NR>1{print $2}' "$B" | sort -u > "$TMP/b.pos"
comm -23 "$TMP/a.pos" "$TMP/b.pos" | awk '{print $1"\tOLD-only"}'  > "$TMP/tag"
comm -13 "$TMP/a.pos" "$TMP/b.pos" | awk '{print $1"\tNEW-only"}' >> "$TMP/tag"
awk '{print "chrX\t"$1}' "$TMP/tag" | sort -k2,2n > "$TMP/rgn"

# 2. whole-chr call at those positions: POS QUAL nALT wmaxDP wtotDP
bcftools query -R "$TMP/rgn" -f '%POS\t%QUAL\t%ALT[\t%AD]\n' "$BCF" 2>/dev/null \
 | awk 'BEGIN{OFS="\t"}{ pos=$1;q=$2; na=($3=="."?0:split($3,aa,","))
        mx=0;tot=0; for(i=4;i<=NF;i++){ if($i=="."||$i==".,.")continue
          n=split($i,ad,","); d=0; for(j=1;j<=n;j++)d+=ad[j]+0; tot+=d; if(d>mx)mx=d }
        print pos,q,na,mx,tot }' > "$TMP/wchr"

# 3. tiled-table per-sample depth (reads-identical check for NEW-only): POS tmax ttot
awk 'FNR>1{ mx=0;tot=0; ns=(NF-2)/2; for(k=1;k<=ns;k++){ d=$(2*k+1)+$(2*k+2); tot+=d; if(d>mx)mx=d }
            print $2"\t"mx"\t"tot }' "$B" > "$TMP/tiled"

# 4. accounting
awk -v TILE="$TILE" -v NC="$NEARCAP" -v T="$TMP/tag" -v L="$TMP/tiled" -v W="$TMP/wchr" '
  FILENAME==T { tag[$1]=$2; next }
  FILENAME==L { tmx[$1]=$2; ttot[$1]=$3; next }
  FILENAME==W { inw[$1]=1; wq[$1]=$2; wn[$1]=$3; wmx[$1]=$4; wtot[$1]=$5; next }
  END{
    OFS="\t"
    print "POS","caller","dL","dR","wQUAL","wnALT","wmaxDP","reads_same","VERDICT","note"
    for(p in tag){
      k=int((p-1)/TILE); cs=k*TILE+1; ce=(k+1)*TILE; dL=p-cs; dR=ce-p
      c=tag[p]; is=(p in inw)
      q=is?wq[p]:"-"; na=is?wn[p]:0; mx=is?wmx[p]:"-"
      near=(is && wmx[p]+0>=NC)
      same="-"; if(c=="NEW-only" && (p in tmx)) same=(is && wmx[p]==tmx[p] && wtot[p]==ttot[p])?"YES":"no"

      if(c=="NEW-only"){
        if(!is){ v="NO-CALL"; nt="whole-chr emitted no variant record here" }
        else if(na>=2){ v="EXPLAINED-3rd"; nt="whole-chr called "na" ALTs -> triallelic -> -m2-M2 dropped it; tiled biallelic/kept" }
        else if(near){ v="EXPLAINED-QUAL"; nt="capped sample (wmaxDP="mx"), biallelic, whole-chr QUAL="q"<=59" }
        else { v=(same=="YES")?"RESIDUAL":"RESIDUAL?"; nt="NO capped sample (wmaxDP="mx"); reads_same="same"; whole-chr QUAL="q"<=59 -> window/QUAL, not subsampling" }
      } else {  # OLD-only: whole-chr kept (QUAL>59); tiled dropped; tiled BCF gone
        if(near){ v="EXPLAINED-cap"; nt="capped sample (wmaxDP="mx") -> subsampling in play in tiled (3rd-allele or QUAL); exact needs tiled BCF" }
        else    { v="RESIDUAL"; nt="NO capped sample (wmaxDP="mx") -> tiled saw identical reads yet dropped it -> window/QUAL" }
      }
      printf "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", p,c,dL,dR,q,na,mx,same,v,nt
    }
  }
' "$TMP/tag" "$TMP/tiled" "$TMP/wchr" > "$TMP/rows"

# header + body sorted by verdict then position, aligned
{ head -1 "$TMP/rows"; tail -n +2 "$TMP/rows" | sort -t$'\t' -k9,9 -k1,1n; } | column -t -s$'\t'

echo
echo "tally by verdict:"
tail -n +2 "$TMP/rows" | cut -f9 | sort | uniq -c | sort -rn
