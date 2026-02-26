#!/bin/bash
conda activate pixyconda

whatshap phase\
	-o /phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04_GQ30DP5_renamed_WHphased.b.vcf.gz\
	--reference Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta\
	/phasing_accuracy/snp_calls/ERR15664662_primary_sorted_ES04_GQ30DP5_renamed.vcf.gz\
	/phasing_accuracy/mapped_reads/ERR15664662_primary_sorted_ES04.RG.bam
