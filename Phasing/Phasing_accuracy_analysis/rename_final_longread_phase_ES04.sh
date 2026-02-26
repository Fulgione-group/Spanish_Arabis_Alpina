#!/bin/bash

conda activate env_bcftools

cd /phasing_accuracy/snp_calls/


bcftools reheader \
  -s <(printf "ES04_hifi\t4\n") \
  -o ES04.long.renamed.isec.vcf.gz \
  ES04.long.isec.vcf.gz

tabix -f -p vcf ES04.long.renamed.isec.vcf.gz




