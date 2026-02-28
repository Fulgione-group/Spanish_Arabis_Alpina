#!/bin/bash
java -jar CalcDiv_DownSample_Allsites_final.jar\
	--vcf GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.b.ParaMask.SNPable.cantabria_france.vcf.gz\
	--cutoff 0.1\
	--binsize 10000\
	--haploidize\
	--popfile unrelated_individuals_new_cantabria_france.popfile
