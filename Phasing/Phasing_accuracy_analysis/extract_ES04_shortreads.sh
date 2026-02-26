#!/bin/bash

conda activate env_bcftools

cd /phasing_accuracy/snp_calls/

bcftools view -s '4' \
  --threads 10\
  GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.intersect.anc_flip_flag.b.vcf.gz \
| bcftools view \
  -v snps -m2 -M2 \
  -g het \
  --threads 10\
  -Oz -o phased_short_reads_ES04.snps.hets.vcf.gz

tabix -f -p vcf phased_short_reads_ES04.snps.hets.vcf.gz

