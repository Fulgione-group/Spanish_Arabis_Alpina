#!/bin/bash
java -jar CreateInput_3PCLR_bugfixed.jar\
	--haploidize\
	--chr $chr\
	--cutoff 0.1\
	--listofpop $Group1,$Group2,$Reference\
	--map geneticMap_${chr}_3pclr_10kb_mindist_rrates.txt\
	--popfile $popfile\
	--vcf Cantabrian_${Group1}_${Group2}_${Reference}_10Permissing_norm_polarized_commonsites.b.vcf.gz
