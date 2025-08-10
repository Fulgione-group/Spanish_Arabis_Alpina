#!/bin/bash
maxjob=${maxjob}
relate='${PATH_TO_RELATE}/relate_v1.2.2_x86_64_static/bin/Relate'


cd $PATH_TO_WD/Relate_no_genes_0kb/bootstraps

anc_mut () {
	local b="$1"
        local c="$2"
	$relate\
		--mode All\
	        -m 7e-9\
	        -N 200000\
	        --haps boot${b}_chr${c}.haps\
	        --sample CAN_FR.sample\
	        --map boot${b}_chr${c}.map\
	        --dist boot${b}_chr${c}.dist\
	        --output relate_GAN_FR_boot${b}_chr${c}
}


for i in {1..100}
do
	for j in {1..8}
	do
	        anc_mut $i $j &
	        background=( $(jobs -p) )
	        if (( ${#background[@]} == $maxjob )); then
	                        wait -n
	        fi
	done
done
wait
