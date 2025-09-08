#!/bin/bash

java -jar Polarize_VCF_from_Fasta_v2.jar \
 --vcf Cantabrian_GroupX_10Permissing.b.vcf.gz\
 --fasta all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_chrall_sorted_mlines_noGaps_uppercases.fasta\
 --stdout | bgzip -c > Cantabrian_GroupX_10Permissing_polarized.b.vcf.gz

tabix Cantabrian_GroupX_10Permissing_polarized.b.vcf.gz
