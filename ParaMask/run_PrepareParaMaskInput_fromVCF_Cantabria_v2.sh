#!/bin/bash
java -jar PrepareParaMaskInput_fromVCF.jar\
        --vcf GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.Cantabrian4Paramask.Clean.30PerMissing.snps.biallelic.b.vcf.gz\
        --missingness 0.3

