#!/bin/bash
plink --bfile CantabrianParaOutSNPable10PercentMissingness\
 --extract CantabrianParaOutSNPable10PercentMissingness_indep.prune.in\
 --make-bed\
 --keep-allele-order\
 --allow-extra-chr\
 --chr 1-8\
 --out CantabrianParaOutSNPable10PercentMissingness_pruned

