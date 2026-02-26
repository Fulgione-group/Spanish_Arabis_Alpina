#!/bin/bash

conda activate env_bcftools

cd /phasing_accuracy/mapped_reads/


samtools addreplacerg \
  -r 'ID:ES04_hifi\tSM:ES04_hifi\tPL:PACBIO\tLB:hifi\tPU:unit1' \
  -o ERR15664662_primary_sorted_ES04.RG.bam \
  ERR15664662_primary_sorted_ES04.bam

samtools index ERR15664662_primary_sorted_ES04.RG.bam

cd /netscratch/dep_coupland/grp_fulgione/bastiaan/data/phasing_accuracy/snp_calls/
bcftools reheader \
  -s <(printf "ERR15664662_primary_sorted_ES04.bam\tES04_hifi\n") \
  -o ERR15664662_primary_sorted_ES04_Q30DP5_renamed.vcf.gz \
  ERR15664662_primary_sorted_ES04_Q30DP5.b.vcf.gz

tabix -f -p vcf ERR15664662_primary_sorted_ES04_Q30DP5_renamed.vcf.gz
