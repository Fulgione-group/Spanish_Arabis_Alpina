# Pipeline to run PCoA using PLINK v1.9 MDS and the R package MASS v7.3-53.1

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

## Run MDS algorithm and plot
- Perform PCoA using **classical multidimensional scaling (`cmdscale`)** from the **MASS v7.3-53.1** package in R:  
  `MDS_Cantabria.R`

- Perform PCoA for Cantabria and France, preprocessing was done the same as for Cantabrians only:
  `MDS_CAN_FR.R` 
