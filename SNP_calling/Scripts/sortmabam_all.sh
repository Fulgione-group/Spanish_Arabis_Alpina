#!/bin/bash

dir=$1
echo "dir"
samples=$2
maxjob=$(expr $3 - 1 )
cd $dir

gatk2='PATH_TO_GATK/gatk-4.2.0.0/gatk'

sortmabam () {
	local j=$1
	Out=$(echo $j | rev | cut -d '_' -f2- | rev)
	$gatk2 --java-options "-XX:ConcGCThreads=2 -XX:ActiveProcessorCount=4" SortSam \
    		--INPUT $j \
    		--OUTPUT ${Out}_mergebamalignment_sorted.bam \
    		--SORT_ORDER coordinate \
    		--TMP_DIR /scratch
}


for i in $(eval ls ${samples}_mergebamalignment.bam)
do
	sortmabam $i & 
	background=( $(jobs -p) )
	if (( ${#background[@]} == $maxjob )); then
		wait -n
	fi
done
wait

