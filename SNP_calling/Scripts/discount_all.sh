#!/bin/bash/

dir=$1
echo "dir"
samples=$2
maxjob=$(expr $3 - 1 )
cd $dir
echo $maxjob
gatk2='PATH_TO_GATK/gatk-4.2.0.0/gatk'


discount_adapters () {
	local j="$1"
	$gatk2 --java-options "-XX:ConcGCThreads=1 -XX:ActiveProcessorCount=${maxjob}" SamToFastq \
        	-I $j \
        	--FASTQ $fastq \
        	--CLIPPING_ATTRIBUTE XT \
        	--CLIPPING_ACTION 2 \
       		--INTERLEAVE true \
        	--NON_PF true
	gzip $fastq
}


for i in $(eval ls ${samples}_adaptermarked.bam)
do
	fastq=$(echo $i | rev | cut -d '_' -f2- | rev)_samtofastq_interleaved.fastq
	discount_adapters $i &
	background=( $(jobs -p) )
      	if (( ${#background[@]} == $maxjob )); then
               		wait -n
	fi
done
wait
