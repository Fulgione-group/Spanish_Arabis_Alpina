#!/bin/bash

conda activate env_bcftools

bedtools intersect -v\
 -a GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.intersect.b.vcf.gz\
 -b ParaMask_SNPable_combined_masks/Arabis_alpina_mpipz_v5.1_annotation.genes.chrall.ZeroBasedStart.bed\
 -header\
 -wa | bgzip -c \
 > GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.intersect.nogenes0kb.b.vcf.gz

