#f!/bin/bash
bcftools view GATK4.2_Allsamples_chrall.filteredQ30LD5UD100K.final.b.vcf.gz\
	-S Cantabrian_GroupX_IDsonly.txt\
	--force-samples\
	| bcftools view -i 'F_MISSING<0.1'\
	| bgzip -c > Cantabrian_GroupX_10Permissing.b.vcf.gz
wait

tabix Cantabrian_GroupX_10Permissing.b.vcf.gz
