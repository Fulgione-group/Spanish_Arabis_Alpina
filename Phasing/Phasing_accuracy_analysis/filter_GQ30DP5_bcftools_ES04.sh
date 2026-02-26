#!/bin/bash

conda activate env_bcftools

cd /phasing_accuracy/mapped_reads/


bcftools view -v snps -m2 -M2 /phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04.b.vcf.gz \
| bcftools filter -i 'QUAL>=30 && DP>=5' \
-Oz -o /phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04_GQ30DP5.b.vcf.gz

tabix -p vcf /phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04_GQ30DP5.b.vcf.gz


