#!/bin/bash
bcftools merge\
 Cantabrian_GroupA_10Permissing_polarized.b.vcf.gz\
 Cantabrian_GroupB_10Permissing_polarized.b.vcf.gz\
 Ref_10Permissing_polarized.b.vcf.gz\
 -Oz\
 -o Cantabrian_GroupA_GroupB_Ref_10Permissing_norm_polarized.b.vcf.gz

