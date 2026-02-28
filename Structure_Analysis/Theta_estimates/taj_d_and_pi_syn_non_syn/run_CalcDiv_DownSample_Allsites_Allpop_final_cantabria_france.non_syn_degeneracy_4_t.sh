#!/bin/bash
java -jar CalcDiv_DownSample_Allsites_final.jar\
	--vcf GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.b.ParaMask.SNPable.cantabria_france.non_syn_degeneracy_4_t.vcf.gz\
	--cutoff 0.1\
	--binsize 10000\
	--haploidize\
	--popfile /netscratch/dep_coupland/grp_fulgione/mehak/metadata/spanish_project_cantabria_france/unrelated_individuals_new_cantabria_france.popfile
