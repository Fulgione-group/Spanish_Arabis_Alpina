#!/bin/bash
conda activate pixyconda

cd /phasing_accuracy/snp_calls

whatshap compare \
  --tsv-pairwise ES04_long_vs_short.tsv \
  ES04.long.renamed.isec.vcf.gz \
  ES04.short.isec.vcf.gz

