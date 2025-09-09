# 3PCLR analysis

java source files and bytecode are in separate directories  "java_src" and "java_bytecode"

## Prepare input

Extract single subject an reference populations with maximum 10% missing genotypes: "Extract_Cantabrian_GroupX.sh"
Polarize each of the subject and reference populations with outgroup alleles: "polarize_3pclr_groupX.sh"
Extract positions for each subject and reference population: "awk '{print $1, $2 - 1, $2}' GroupX_polarized.b.vcf.gz > GroupX_ZeroBasedStart.bed" 
Merge two subjects and reference population: "merge_4_3PCLR_GroupA_GroupB_Ref.sh"
Filter the merged VCF file to retain sites where all subjects and refs have maximum of 10% missing genotypes: "filter_sites4Driftrates_GroupA_GroupB_Ref.sh"

Create input from a VCF file per chromosome using a genetic map, missingness cutoff of 10%, two subject populations (GroupA, GroupB) and one reference population (Ref), and sampling one haplotype per genotype per site: "run_CreatInput_3PCLR_chrX_GroupA_GroupB_Ref.sh"

## Extra steps to calculate drift rates

Apply ParaMasks for Cantabria and France, apply SNPable mask, and 2kbp Up and downstream of genes: "filter_bedtools_Can4Drift_SNPable_CanParaOut_FranParaOut_noGenes_2kb_GroupA_GroupB_Ref.sh"
Make 3d SFS with projecting the populations down by 10% and sampling only one haplotype per genotype per site: "run_Downsample_4Driftrates_GroupA_GroupB_Ref.sh"
Make input table for drift rate calculation script, Nhap correspond to the number of haplotypes: "awk 'BEGIN{FS=OFS="\t"}{print $1, Nhap_GroupA, $2, Nhap_GroupB, $3,Nhap_Ref, $4}' GroupA_GroupB_Ref.SFS > GroupA_GroupB_Ref_SFS.tab"
Calculate F3 stats: "Rscript CalcDriftsF3.R GroupA_GroupB_Ref_SFS.tab GroupA_GroupB_Ref_SFS.F3"

## 3PCLR

run 3PCLR analysis per chromosome with minimum distance of 25 bp , 100 snps per window, : "run_3pclr_G1_G2_REF_chrX.sh"
