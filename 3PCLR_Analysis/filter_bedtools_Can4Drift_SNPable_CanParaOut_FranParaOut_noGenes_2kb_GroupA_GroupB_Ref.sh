#!/bin/bash

bedtools intersect -v\
 -a Cantabrian_GroupA_GroupB_Ref_10Permissing_polarized_commonsites.b.vcf.gz\
 -b ParaMask_Cantabria_and_French_30PerMiss_v3_EM2.7.1_EMresults.chrall.finalClass.multicopy.ZeroBasedStart.SNPable.3col.sorted.merged.bed.noGenes2kbpUpandDown.sorted.merged.bed\
 -header\
 -wa | bgzip -c \
 > Cantabrian_GroupA_GroupB_Ref_10Permissing_polarized_commonsites_CanParaOut_FranParaOut_noGenes2kb.b.vcf.gz

wait
tabix Cantabrian_GroupA_GroupB_Ref_10Permissing_polarized_commonsites_CanParaOut_FranParaOut_noGenes2kb.b.vcf.gz
