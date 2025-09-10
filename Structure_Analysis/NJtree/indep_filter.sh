#!/bin/bash
plink --bfile CantabrianParaOutSNPable10PercentMissingness\
 --indep-pairwise 50 10 0.1\
 --keep-allele-order\
 --allow-extra-chr\
 --out CantabrianParaOutSNPable10PercentMissingness_indep
