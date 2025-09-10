#!/bin/bash
plink --vcf GATK4.2_All_chrall.filteredQ30LD5UD100K.final.Cantabrian4Paramask.Clean.10PerMissing.snps.biallelic.unrelatedNew.ParaOut.SNPable.NoSingleton.ID.b.vcf.gz\
 --biallelic-only --snps-only\
 --make-bed\
 --mac 2\
 --keep-allele-order\
 --allow-extra-chr\
 --geno 0.1\
 --out CantabrianParaOut10PercentMissingnessNoSingletonsUnrelated

