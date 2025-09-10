# Pipeline to make Neighbour-Joining tree using plink and ape package
Prior to the analyis only sites with less than 10%  missing genotypes were retained. Multicopy regions inferred by ParaMask and unmappable regions inferred by SNPable.

## Pre-process
Covert clean vcf to bed file: "make_bed.sh"
Make indep filter for a sliding window of 50 bp, stepsize of 10 bp and maximum linkage of 0.1 : "indep_filter.sh"
Prune SNPs from indep filter: "prune_ld.sh"
create 1-IBS distance matrix: "make_distance_matrix.sh"

## Create NJ tree and plot

run NJ algorithm from ape and plot: "NJ_Cantabria.R"

