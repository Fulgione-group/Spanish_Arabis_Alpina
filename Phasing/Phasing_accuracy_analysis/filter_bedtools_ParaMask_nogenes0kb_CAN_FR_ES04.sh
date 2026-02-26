#!/bin/bash

conda activate env_bcftools

bedtools intersect -v\
 -a /phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased.het.b.vcf.gz\
 -b ParaMask_SNPable_combined_masks/ParaMask_Cantabria_and_French_30PerMiss_v3_EM2.7.1_EMresults.chrall.finalClass.multicopy.ZeroBasedStart.SNPable.3col.sorted.merged.bed\
 -header\
 -wa | bgzip -c \
 > /phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased_ParaOut_CANFR_snpable.het.b.vcf.gz

