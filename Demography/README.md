#Demographic Analysis using relate


##Point estimates

###infer genealogies using relate:
for demography: "anc_mut_chr_CAN_FR_no_genes_0kb.sh"
for genic regions: "anc_mut_chr_CAN_FR_with_genes.sh"

###update coalescence rates and estimate popsizes
for demography: "EstimatePopSize_CAN_FR_no_genes_0kb.sh"
for genic regions: "EstimatePopSize_CAN_FR_with_genes.sh"



##bootstrapping

###create bootstraps
divide the genome in 5 Mbp segments and resample to full genomic size: "relate_popsize_CAN_FR_boot_make_chr.R"

###infer genealogies using relate
"anc_mut_boot_all_0kbp.sh"

###update coalecence rates and estimate popsizes
"Estimate_popsize_boot.sh"



##Downstream analysis

###convert to treeseqeuence format 
"Convert_relateWithGenes_to_treeseq_chr_final.sh"

###Extract genealogy at FRL1 non-synonymous mutations
"Extract_FRIL1_mut_12619412.sh"
"Extract_FRIL1_mut_12619091.sh"
"Extract_FRIL1_mut_12619062.sh"
"Extract_FRIL1_mut_12619056.sh"
"Extract_FRIL1_mut_12618944.sh"
"Extract_FRIL1_mut_12618926.sh"
