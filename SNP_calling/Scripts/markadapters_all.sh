#!/bin/bash/

dir=$1
echo "dir"
samples=$2
maxjob=$3
cd $dir

gatk2='PATH_TO_GATK/gatk-4.2.0.0/gatk'


markadapters () {
	local j="$1"
	gunzip -k $j
	uj=$(echo $j | rev | cut -d '.' -f2- | rev)
	OUT=$(echo ${uj} | rev | cut -d '_' -f2- | rev)
	$gatk2 --java-options "-XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:ActiveProcessorCount=1" MarkIlluminaAdapters \
        	--INPUT ${uj} \
        	--OUTPUT ${OUT}_adaptermarked.bam \
        	--METRICS ${OUT}_metrics.txt \
		--TMP_DIR /scratch\
		--FIVE_PRIME_ADAPTER AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
		--THREE_PRIME_ADAPTER AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGA
	rm ${uj}
}

for i in $(eval ls ${samples}_fastqtosam.sam.gz)
do
	markadapters $i &
	background=( $(jobs -p) )
	if (( ${#background[@]} == $maxjob )); then
        	wait -n
	fi
done
