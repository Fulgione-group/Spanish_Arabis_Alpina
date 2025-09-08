#!/bin/bash

bedtools intersect\
 -a Cantabrian_GroupA_GroupB_Ref_10Permissing_norm_polarized.b.vcf.gz\
 -b common_GroupA_GroupB_Ref_sites_pos_ZeroBasedStart.bed\
 -header\
 -wa | bgzip -c > Cantabrian_GroupA_GroupB_Ref_10Permissing_norm_polarized_commonsites.b.vcf.gz

tabix Cantabrian_GroupA_GroupB_Ref_10Permissing_norm_polarized_commonsites.b.vcf.gz

