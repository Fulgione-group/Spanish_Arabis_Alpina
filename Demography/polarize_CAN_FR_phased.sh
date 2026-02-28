#!/bin/bash
java -jar Polarize_VCF_from_Fasta.jar \
 --vcf GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.b.vcf.gz\
 --fasta all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_chrall_sorted_mlines_noGaps_uppercases.fasta\
 --stdout | bgzip -c > GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.pall.CAN.FR.prephased.20PercentMissing.ParaOut.SNPable.phased.chrall.allsites.polarized.b.vcf.gz

