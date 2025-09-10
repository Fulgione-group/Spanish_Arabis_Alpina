#!/bin/bash

java -jar Calc_Div_F_projectDown.jar\
	--vcf GATK4.2_All_chrall.filteredQ30LD5UD100K.final.Cantabrian4Paramask.Clean.10PerMissing.snps.biallelic.ParOut.SNPable.ID.b.vcf.gz\
	--cutoff 0.1\
	--popfile unrelated_individuals_new_4Fis.popfile
