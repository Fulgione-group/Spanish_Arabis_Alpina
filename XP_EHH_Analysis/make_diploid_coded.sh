#!/bin/bash
zcat GROUP1.haploidized.vcf.gz | \
awk 'BEGIN{OFS="\t"}
/^#/ {print; next}
{
  for(i=10;i<=NF;i++){
    if($i=="0") $i="0|0";
    else if($i=="1") $i="1|1";
    else $i=".|.";
  }
  print
}' | bgzip > GROUP1.haploidized.diploidcoded.vcf.gz

tabix -p vcf GROUP1.haploidized.diploidcoded.vcf.gz

