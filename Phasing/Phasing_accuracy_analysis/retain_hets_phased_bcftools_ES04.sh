#!/bin/bash

conda activate env_bcftools

cd /phasing_accuracy/snp_calls/


bcftools view -v snps -m2 -M2 -g het ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased.b.vcf.gz \
| bcftools filter -i 'QUAL>=30 && DP>=5' \
-Oz -o ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased.het.b.vcf.gz

tabix -p vcf ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased.het.b.vcf.gz


