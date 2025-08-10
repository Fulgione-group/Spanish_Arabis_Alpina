#!/bin/bash
maxjob=100
popsize='${PATH_TO_RELATE}/relate_v1.2.2_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh'


cd $PATH_TO_WD/Relate_no_genes/bootstraps

estimate_ps () {
	local b="$1"
	$popsize\
		-i relate_GAN_FR_boot${b} \
		--first_chr 1 \
		--last_chr 8 \
		-m 7e-9 \
		--years_per_gen 1.5 \
		--poplabels CAN_FR.poplabels \
		--threads 20 \
		--num-iter 10 \
		-o relate_popsize_CAN_FR_boot${b}
}


for i in {1..100}
do
        estimate_ps $i &
        background=( $(jobs -p) )
        if (( ${#background[@]} == $maxjob )); then
                        wait -n
        fi
done
wait
