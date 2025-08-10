#!/bin/bash
cd $PATH_TO_WD/Relate_with_genes

$PATH_TO_RELATE/relate_v1.2.2_x86_64_static/bin/Relate\
      	--mode All\
      	-m 7e-9\
	-N 200000\
	--haps GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.intersect.chr1.snps.haps\
        --sample CAN_FR.sample\
      	--map /geneticMap_${chr}_4SpanishProject_relate.map\
	--dist GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.intersect.${chr}.dist\
        --output relate_GAN_FR_${chr}

