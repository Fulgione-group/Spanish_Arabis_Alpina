#!/bin/bash
cd $PATH_TO_WD/Relate_with_genes

$PATH_TO_RELATE/relate_v1.2.2_x86_64_static/scripts/TreeView/TreeViewMutation.sh \
                 --haps GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.intersect.chr8.snps.haps \
                 --sample CAN_FR.sample \
                 --anc relate_popsize_CAN_FR_chr8.anc.gz \
                 --mut relate_popsize_CAN_FR_chr8.mut.gz \
                 --poplabels CAN_FR_FRIL1.poplabels \
                 --bp_of_interest 12618944  \
                 --years_per_gen 1.5 \
                 -o FRIL1_mutation_12618944
