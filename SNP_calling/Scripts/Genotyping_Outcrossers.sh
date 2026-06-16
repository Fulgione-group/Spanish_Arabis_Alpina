#!/bin/bash
set -euo pipefail

ref="Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta"
positions="positions.txt"

# positions.txt: chrom<TAB>pos
awk '{print $1"\t"$2-1"\t"$2}' "$positions" > positions.bed

bcftools mpileup \
    -f "$ref" \
    -q 30 \
    -R positions.bed \
    -Ou \
    *_gatk42_bwamem.sorted.bam | \
bcftools call \
    -mv \
    -Oz \
    -o mpileup_positions.vcf.gz

bcftools index mpileup_positions.vcf.gz
