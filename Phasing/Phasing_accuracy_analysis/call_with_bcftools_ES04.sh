#!/bin/bash

conda activate env_bcftools

cd /phasing_accuracy/mapped_reads/

bcftools mpileup \
  -Ou \
  -f Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta \
  ERR15664662_primary_sorted_ES04.bam \
| bcftools call \
  -mv \
  -Oz \
  -o /phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04.b.vcf.gz
