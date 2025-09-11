# Admixture Pipeline

Prior to the admixture analysis, only unrelated samples and sites with less than 10% missing genotypes were retained.  
Multicopy regions inferred by ParaMask and unmappable regions inferred by SNPable, as well as singletons, were filtered out.

## Pre-process
- Convert clean VCF to BED file:  
  `make_bed.sh`

- Make LD independence filter with a sliding window of 50 SNPs, step size of 10 SNPs, and maximum linkage threshold of 0.1:  
  `indep_filter.sh`

- Prune SNPs from independence filter:  
  `prune_ld.sh`

## Run ADMIXTURE
- Admixture analysis with up to 20 ancestry groups (repeated 10 times) using **ADMIXTURE v1.3.0** (Alexander et al., 2009).  
  Cross-validation (`--cv`) was used to evaluate the best-supported number of ancestry groups (*K*):  
  `ADMIXTURE_nostderr_1.sh`

## Run pong
- Find the best alignment of ancestry groups between runs using **pong** (Behr et al., 2016):  
  `pong.sh`
