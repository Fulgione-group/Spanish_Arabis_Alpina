# Admixture Pipeline

Prior to the admixture analysis only unrelated samples, and only sites with less than 10%  missing genotypes were retained. Multicopy regions inferred by ParaMask and unmappable regions inferred by SNPable as well as singletons were filtered out.

## Pre-process
Covert clean vcf to bed file: "make_bed.sh"
Make indep filter for a sliding window of 50 bp, stepsize of 10 bp and maximum linkage of 0.1 : "indep_filter.sh"
Prune SNPs from indep filter: "prune_ld.sh"

## run admixture
Admixture with up to 20 number of ancestry groups (this was repeated 10 times): "ADMIXTURE_nostderr_1.sh"

## run Pong 
Find best alignment of ancestry groups between runs: "pong.sh"
