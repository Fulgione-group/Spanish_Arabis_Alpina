#!/bin/bash
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8
do
  bcftools view -r $chr GROUP1.haploidized.diploidcoded.vcf.gz -Oz -o GROUP1.${chr}.vcf.gz
  bcftools view -r $chr GROUP2.haploidized.diploidcoded.vcf.gz -Oz -o GROUP2.${chr}.vcf.gz
done
