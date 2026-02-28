# Calculation of Watterson’s theta, nucleotide diversity (π), and Tajima’s D genome-wide

Java source files and bytecode for the custom scripts are located in the subdirectories.

Prior to the analysis, only unrelated samples and sites with less than 10% missing genotypes were retained.  
Multicopy regions inferred by ParaMask, unmappable regions inferred by SNPable, as well as singletons, were filtered out.

Samples per population were projected down to 90% of the total sample size (to account for missing data) from a VCF that includes invariant sites.  
Including invariant sites ensures preservation of SNP density and correct denominators in the calculations.  
Genotypes and populations are specified in the `popfile`.

The program outputs both **bin-wise estimates** (per genomic window) and **genome-wide estimates** of Watterson’s theta, nucleotide diversity (π), and Tajima’s D, scripts for calculations on synonymous and non-synonymous sites **taj_d_and_pi_syn_non_syn**.

- Run calculation:  
  `run_CalcDiv_DownSample_Allsites_Allpop_final.sh`


