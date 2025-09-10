#!/bin/bash
java -jar CalcDiv_DownSample_Allsites_final.jar\
	--vcf GATK4.2_All_chrall.filteredQ30LD5UD100K.final.Cantabrian.Clean.10PerMissing.Allsites.biallelic.updated.unrelatedNew.b.vcf.gz\
	--cutoff 0.1\
	--binsize 10000\
	--haploidize\
	--popfile unrelated_individuals_new.popfile
