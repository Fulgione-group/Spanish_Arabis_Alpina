# Calculation of Theta watterson and pi, and TajD genome-wide

java src and bytcode for the custom scripts are in the subdirectories

Prior to the admixture analysis only unrelated samples, and only sites with less than 10%  missing genotypes were retained. Multicopy regions inferred by ParaMask and unmappable regions inferred by SNPable as well as singletons were filtered out.
Project down samples per population (to 1-cutoff = 90%) from a VCF that includes invariant sites (to preserve SNP densities and ensure correct denumerators), genotypes and populations are specified in the popfile: "run_CalcDiv_DownSample_Allsites_Allpop_final.sh"

