#!/bin/bash
cd $PATH_TO_WD/Relate_with_genes

$PATH_TO_RELATE_LIB/relate_lib/bin/Convert\
              --mode ConvertToTreeSequence \
	      --compress \
              --anc relate_popsize_CAN_FR_${chr}.anc.gz \
              --mut relate_popsize_CAN_FR_${chr}.mut.gz \
              -o relate_popsize_CAN_FR_msp_final_${chr}
