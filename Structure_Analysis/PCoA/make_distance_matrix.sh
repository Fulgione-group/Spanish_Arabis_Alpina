#!/bin/bash
source /netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/anaconda3/etc/profile.d/conda.sh
conda activate baconda



plink --bfile CantabrianParaOutSNPable10PercentMissingness_pruned\
 --out CantabrianParaOutSNPable10PercentMissingness_pruned_clustering\
 --distance '1-ibs'

