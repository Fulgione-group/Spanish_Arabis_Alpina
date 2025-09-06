#!/bin/bash
Rscript --vanilla ParaMask_EM_v0.2.7.1.R\
      	--het GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.Cantabrian4Paramask.Clean.30PerMissing.snps.biallelic.b.vcf.het.stat.txt\
        --missingness 0.3\
	--outdir ./\
	--ID Cantabria_30PerMiss_v3_EM2.7.1
