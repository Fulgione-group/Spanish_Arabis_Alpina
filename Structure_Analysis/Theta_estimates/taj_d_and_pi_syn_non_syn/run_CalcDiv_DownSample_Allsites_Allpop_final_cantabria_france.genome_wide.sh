#!/bin/bash
java -jar CalcDiv_DownSample_Allsites_final.jar\
	--vcf /netscratch/dep_coupland/grp_fulgione/mehak/data/masks_filtered/paramask_snpable_filetered/GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.b.ParaMask.SNPable.cantabria_france.vcf.gz\
	--cutoff 0.1\
	--binsize 10000\
	--haploidize\
	--popfile /netscratch/dep_coupland/grp_fulgione/mehak/metadata/spanish_project_cantabria_france/unrelated_individuals_new_cantabria_france.popfile
