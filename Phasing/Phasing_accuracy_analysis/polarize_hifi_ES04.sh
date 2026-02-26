#!/bin/bash
source /netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/anaconda3/etc/profile.d/conda.sh
conda activate env_bcftools


java -jar Polarize_VCF_from_Fasta_v2.jar \
 --vcf /phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased_ParaOut_CANFR_snpable.het.b.vcf.gz\
 --fasta all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_chrall_sorted_mlines_noGaps_uppercases.fasta\
 --stdout | bgzip -c > \
/phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased_ParaOut_CANFR_snpable.het.polarized.b.vcf.gz

tabix /phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased_ParaOut_CANFR_snpable.het.polarized.b.vcf.gz
