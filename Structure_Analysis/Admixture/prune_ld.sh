#!/bin/bash
source /netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/anaconda3/etc/profile.d/conda.sh
conda activate baconda

plink --bfile CantabrianParaOut10PercentMissingnessNoSingletonsUnrelated\
 --extract CantabrianParaOut10PercentMissingnessNoSingletonsUnrelated_indep.prune.in\
 --make-bed\
 --keep-allele-order\
 --allow-extra-chr\
 --chr 1-8\
 --out CantabrianParaOut10PercentMissingnessNoSingletonsUnrelated_pruned

