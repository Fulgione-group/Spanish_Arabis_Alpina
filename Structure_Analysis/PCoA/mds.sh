#!/bin/bash
plink --allow-extra-chr\
 --bfile\
 CantabrianParaOutSNPable10PercentMissingness_pruned\
 --cluster\
 --mds-plot 4 eigvals
