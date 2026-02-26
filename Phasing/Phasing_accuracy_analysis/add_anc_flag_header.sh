#!/bin/bash
conda activate env_bcftools

cd /phasing_accuracy/snp_calls


bcftools annotate -h add_filter_ancflip.hdr \
  --threads 20\
  -Oz -o GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.intersect.anc_flip_flag.b.vcf.gz \
  GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.intersect.b.vcf.gz
