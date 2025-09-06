#!/bin/bash
Rscript --vanilla ParaMask_EM_v0.2.7.1.R\
      	--het GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.French4Paramask.30PerMissing.snps.biallelic.b.vcf.het.stat.txt\
        --missingness 0.3\
	--outdir ./\
	--ID ParaMask_French_clean_EM_2.7.1
