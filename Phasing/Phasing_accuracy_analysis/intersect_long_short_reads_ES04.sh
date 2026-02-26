#!/bin/bash

conda activate env_bcftools

cd /data/phasing_accuracy/snp_calls/


bcftools isec -n=2 -w1 -Oz -o ES04.long.isec.vcf.gz  ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased_ParaOut_CANFR_snpable.het.polarized.onlyphased.anc_flip_flag.b.vcf.gz phased_short_reads_ES04.snps.hets.vcf.gz
bcftools isec -n=2 -w2 -Oz -o ES04.short.isec.vcf.gz ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased_ParaOut_CANFR_snpable.het.polarized.onlyphased.anc_flip_flag.b.vcf.gz phased_short_reads_ES04.snps.hets.vcf.gz
tabix -f -p vcf ES04.long.isec.vcf.gz
tabix -f -p vcf ES04.short.isec.vcf.gz
