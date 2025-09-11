# 3P-CLR analysis

Java source files and bytecode are stored in separate directories: `java_src` and `java_bytecode`.  
The 3P-CLR program is used from Racimo (2016).

## Prepare input

- Extract single subject and reference populations with a maximum of 10% missing genotypes:  
  `Extract_Cantabrian_GroupX.sh`

- Polarize each of the subject and reference populations with outgroup alleles:  
  `polarize_3pclr_groupX.sh`

- Extract positions for each subject and reference population:  
  `awk '{print $1, $2 - 1, $2}' GroupX_polarized.b.vcf.gz > GroupX_ZeroBasedStart.bed`

- Merge two subjects and the reference population:  
  `merge_4_3PCLR_GroupA_GroupB_Ref.sh`

- Filter the merged VCF file to retain sites where all subjects and references have at most 10% missing genotypes:  
  `filter_sites4Driftrates_GroupA_GroupB_Ref.sh`

- Create input from a VCF file per chromosome using a genetic map, missingness cutoff of 10%, two subject populations (GroupA, GroupB), one reference population (Ref), and sampling one haplotype per genotype per site:  
  `run_CreatInput_3PCLR_chrX_GroupA_GroupB_Ref.sh`

## Extra steps to calculate drift rates

- Apply ParaMasks for Cantabria and France, apply the SNPable mask, and exclude 2 kbp up- and downstream of genes:  
  `filter_bedtools_Can4Drift_SNPable_CanParaOut_FranParaOut_noGenes_2kb_GroupA_GroupB_Ref.sh`

- Make 3D SFS by projecting the populations down by 10% and sampling only one haplotype per genotype per site:  
  `run_Downsample_4Driftrates_GroupA_GroupB_Ref.sh`

- Make an input table for the drift rate calculation script (`Nhap` corresponds to the number of haplotypes):  
  `awk 'BEGIN{FS=OFS="\t"}{print $1, Nhap_GroupA, $2, Nhap_GroupB, $3, Nhap_Ref, $4}' GroupA_GroupB_Ref.SFS > GroupA_GroupB_Ref_SFS.tab`

- Calculate F3 statistics:  
  `Rscript CalcDriftsF3.R GroupA_GroupB_Ref_SFS.tab GroupA_GroupB_Ref_SFS.F3`

## 3P-CLR

- Run 3P-CLR analysis per chromosome with a minimum distance of 25 bp and 100 SNPs per window:  
  `run_3pclr_G1_G2_REF_chrX.sh`

## Analysis of 3P-CLR

- Identify overlapping genes in the tails and perform GO analysis of genes in the tails:  
  `3PCLR_final_analysis_V3.R`



