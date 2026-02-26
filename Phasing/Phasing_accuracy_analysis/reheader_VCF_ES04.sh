#!/bin/bash

conda activate env_bcftools


cd /phasing_accuracy/snp_calls/
bcftools reheader \
  -s <(printf "ERR15664662_primary_sorted_ES04.bam\tES03_hifi\n") \
  -o ERR15664662_primary_sorted_ES04_GQ30DP5_renamed.vcf.gz \
  ERR15664662_primary_sorted_ES04_GQ30DP5.b.vcf.gz

tabix -f -p vcf ERR15664662_primary_sorted_ES04_GQ30DP5_renamed.vcf.gz
