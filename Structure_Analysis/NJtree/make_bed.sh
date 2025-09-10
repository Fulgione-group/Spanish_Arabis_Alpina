#!/bin/bash
plink --vcf ../GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.Cantabrian4Paramask.Clean.10PerMissing.snps.biallelic.ParOut.SNPable.ID.b.vcf.gz\
 --biallelic-only --snps-only\
 --make-bed\
 --keep-allele-order\
 --allow-extra-chr\
 --geno 0.1\
 --out CantabrianParaOutSNPable10PercentMissingness

