#!/bin/bash
conda activate env_bcftools


java -jar ~/scripts/java/Polarize_VCF_from_Fasta_v2.jar \
 --vcf 3P_XPCLR_v3/Cantabrian_Group2_10Permissing_norm.b.vcf.gz\
 --fasta all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_chrall_sorted_mlines_noGaps_uppercases.fasta\
 --stdout | bgzip -c > Cantabrian_Group2_10Permissing_norm_polarized.b.vcf.gz

tabix Cantabrian_Group2_10Permissing_norm_polarized.b.vcf.gz
