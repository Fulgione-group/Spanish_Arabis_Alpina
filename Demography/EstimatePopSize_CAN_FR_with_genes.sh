#!/bin/bash
cd $PATH_TO_WD/Relate_with_genes

$PATH_TO_RELATE/relate_v1.2.2_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i relate_GAN_FR \
              --first_chr 1 \
              --last_chr 8 \
              -m 7e-9 \
	      --years_per_gen 1.5 \
	      --poplabels CAN_FR.poplabels \
	      --threads 20 \
	      --num-iter 10\
              -o relate_popsize_CAN_FR
