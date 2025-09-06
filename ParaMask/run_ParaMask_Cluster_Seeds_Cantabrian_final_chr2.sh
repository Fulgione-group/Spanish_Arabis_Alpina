#!/bin/bash
java -jar ParaMask_Cluster_Seeds.jar\
	--cov GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.Cantabrian4Paramask.Clean.30PerMissing.snps.biallelic.b.vcf.het.cov.stat.chr2.txt\
	--het Cantabria_30PerMiss_v3_EM2.7.1_EMresults.chr2.het\
	--covgw  GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.Cantabrian4Paramask.Clean.30PerMissing.snps.biallelic.b.vcf.cov.gw.txt\
	--cutoff  148\
	--range 1,39892743
