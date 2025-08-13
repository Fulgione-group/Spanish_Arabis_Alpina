#!/bin/bash/
module load shapeit4/v4.2.0
cd $PATH_TO_WD
shapeit4.2 --input GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.FR.prephased.20PercentMissing.ParaOut.SNPable.b.vcf.gz\
	--region $chr \
	--use-PS 0.0001 \
	--map geneticMap_${chr}_4SpanishProject_phasing.txt \
	--thread 10 \
	--output GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.${chr}.b.vcf.gz
