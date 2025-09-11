# Pipeline to construct a Neighbor-Joining tree using PLINK v1.9 and the R package ape v5.5

Prior to the analysis, only sites with less than 10% missing genotypes were retained.  
Multicopy regions inferred by ParaMask and unmappable regions inferred by SNPable were filtered out.

## Pre-process
- Convert clean VCF to BED file:  
  `make_bed.sh`

- Make LD independence filter with a sliding window of 50 SNPs, step size of 10 SNPs, and maximum linkage threshold of 0.1:  
  `indep_filter.sh`

- Prune SNPs from independence filter:  
  `prune_ld.sh`

- Create a 1-IBS distance matrix with **PLINK v1.9** (`--distance 1-ibs`):  
  `make_distance_matrix.sh`

## Create NJ tree and plot
- Run the Neighbor-Joining algorithm from the **ape v5.5** package in R and plot the tree:  
  `NJ_Cantabria.R`
