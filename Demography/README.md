# Demographic Analysis using Relate v.1.2.2

## Point estimates

### Infer genealogies with Relate
- For demography:  
  `anc_mut_chr_CAN_FR_no_genes_0kb.sh`  
- For genic regions:  
  `anc_mut_chr_CAN_FR_with_genes.sh`

### Update coalescence rates and estimate population sizes
- For demography:  
  `EstimatePopSize_CAN_FR_no_genes_0kb.sh`  
- For genic regions:  
  `EstimatePopSize_CAN_FR_with_genes.sh`

## Bootstrapping

### Create bootstraps
- Divide the genome into 5 Mbp segments and resample to full genomic size:  
  `relate_popsize_CAN_FR_boot_make_chr.R`

### Infer genealogies with Relate
- `anc_mut_boot_all_0kbp.sh`

### Update coalescence rates and estimate population sizes
- `Estimate_popsize_boot.sh`

## Downstream analysis

### Infer split times and plot results

- `relate_popsize_CAN_FR.R`

### Convert to tree sequence format
- `Convert_relateWithGenes_to_treeseq_chr_final.sh`
  
### Extract genealogies at FRL1 non-synonymous mutations
- `Extract_FRIL1_mut_12619412.sh`  
- `Extract_FRIL1_mut_12619091.sh`  
- `Extract_FRIL1_mut_12619062.sh`  
- `Extract_FRIL1_mut_12619056.sh`  
- `Extract_FRIL1_mut_12618944.sh`  
- `Extract_FRIL1_mut_12618926.sh`
