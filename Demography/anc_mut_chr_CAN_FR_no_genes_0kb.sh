#!/bin/bash
cd $PATH_TO_WD/Relate_no_genes_0kb

$PATH_TO_RELATE/relate_v1.2.2_x86_64_static/bin/Relate\
      	--mode All\
      	-m 7e-9\
	-N 200000\
	--haps GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.intersect.nogenes0kbUpDownStream.${chr}.snps.haps\
        --sample CAN_FR.sample\
      	--map /netscratch/dep_coupland/grp_fulgione/bastiaan/data/genetic_maps/geneticMap_${chr}_4SpanishProject_relate.map\
	--dist GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.intersect.nogenes0kbUpDownStream.${chr}.dist\
        --output relate_CAN_FR_${chr}

