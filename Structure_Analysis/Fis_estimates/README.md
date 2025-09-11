# Calculation of multi-population FIS estimate per bin

To calculate FIS, we used populations with sample sizes of at least 10 to avoid high error rates in allele frequency estimates.  
We calculated Nei’s observed heterozygosity (Ho) and expected heterozygosity (Hs) per site across populations (Nei, 1987).  
Genome-wide and window-based (10,000 SNPs) FIS estimates were then computed as:

FIS = 1 − ΣHo / ΣHs

To account for missing values, we randomly projected the number of genotypes down to 90% of the total sample size per population.

Java source files and bytecode for the custom scripts are in the subdirectories.

Prior to the analysis, only unrelated samples and sites with less than 10% missing genotypes were retained.  
Multicopy regions inferred by ParaMask, unmappable regions inferred by SNPable, as well as singletons, were filtered out.

- Project down samples per population (cutoff = 90%) from a VCF that includes invariant sites  
  (to preserve SNP densities and ensure correct denominators).  
  Genotypes and populations are specified in the `popfile`:  
  `run_Calc_Div_F_projectDown_Cantabria_new.sh`
