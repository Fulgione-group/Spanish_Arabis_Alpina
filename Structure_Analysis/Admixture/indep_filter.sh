#!/bin/bash
plink --bfile CantabrianParaOut10PercentMissingnessNoSingletonsUnrelated\
 --indep-pairwise 50 10 0.1\
 --keep-allele-order\
 --allow-extra-chr\
 --out CantabrianParaOut10PercentMissingnessNoSingletonsUnrelated_indep
