#!/bin/bash
java -jar FilterVCF.jar --vcf Alpina_Genomes.vcf.gz\
	--indels-exclude\
	--biallelic\
	--quality 30\
	--depth 5,100000\
	--coverage GATK4.2_alpina_15_02_2023_chrall.vcf.avcov.txt,-1\
	--min-sample-cov 10\
	--gzip\
	--out Alpina_Genomes.filteredQ30LD5UD100K.vcf

