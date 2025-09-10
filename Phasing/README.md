# Two-step Phasing

## Prephasing
Use mapped BAM files to identify alleles linked on read pairs: "whatshap_GATK4.2_All_SNPS_pall.sh"


## Phasing
Filter out ParaMask, SNPable mask, and sites with more than 20 percent missingness.
Use prephased VCF files to phase and impute per region with shapeit4.
France: "shapeit4_FR_ParaOut_SNPable_20PerMiss_chrX.sh"
Cantabria: "shapeit4_CAN_ParaOut_SNPable_20PerMiss_chrX.sh"
