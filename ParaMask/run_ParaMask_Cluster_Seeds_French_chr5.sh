#!/bin/bash
java -jar ParaMask_Cluster_Seeds.jar\
	--cov GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.French4Paramask.30PerMissing.snps.biallelic.b.vcf.het.cov.stat.chr5.txt\
	--het ParaMask_French_clean_EM_2.7.1_EMresults.chr5.het\
	--covgw  GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.French4Paramask.30PerMissing.snps.biallelic.b.vcf.cov.gw.txt\
	--cutoff  127\
	--range 1,28376109
