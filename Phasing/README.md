# Two-step Phasing

## Prephasing
Use mapped BAM files to identify alleles linked on read pairs with **WhatsHap v1.1**:  
`whatshap_GATK4.2_All_SNPS_pall.sh`

## Phasing
Filter out ParaMask regions, SNPable mask, and sites with more than 20% missingness.  
Use the prephased VCF files to phase and impute per region with **SHAPEIT v4.2**.

- France:  
  `shapeit4_FR_ParaOut_SNPable_20PerMiss_chrX.sh`

- Cantabria:  
  `shapeit4_CAN_ParaOut_SNPable_20PerMiss_chrX.sh`

## Phasing efficiency analysis

- Longread based assesment of phasing accuracy including README in:
`Phasing_accuracy_analysis`
