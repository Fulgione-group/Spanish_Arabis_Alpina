#!/bin/bash
set -euo pipefail

VCF_PREFIX="GROUP1"   # expects GROUP1.chr1.vcf.gz etc.
MAP_PREFIX="geneticMap"
OUT_PREFIX="selscan"

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8
do
    VCF="${VCF_PREFIX}.${chr}.vcf.gz"
    RECMAP="geneticMap_${chr}_3pclr_10kb_mindist_rrates.txt"
    OUT="${OUT_PREFIX}.${chr}.map"

    echo "Making map for ${chr}"

    bcftools query -f '%CHROM\t%POS\n' "$VCF" | \
    awk -v map="$RECMAP" '
    BEGIN{
        OFS="\t"
        while((getline < map) > 0){
            c=$1; pos=$2; morg=$3
            n++
            p[n]=pos
            m[n]=morg
        }
    }
    {
        chr=$1
        pos=$2

        if(pos <= p[1]){
            morg = m[1] * pos / p[1]
        }
        else if(pos >= p[n]){
            morg = m[n]
        }
        else{
            for(i=1; i<n; i++){
                if(pos >= p[i] && pos <= p[i+1]){
                    frac = (pos-p[i])/(p[i+1]-p[i])
                    morg = m[i] + frac*(m[i+1]-m[i])
                    break
                }
            }
        }

        print chr, chr ":" pos, morg*100, pos
    }' > "$OUT"
done
