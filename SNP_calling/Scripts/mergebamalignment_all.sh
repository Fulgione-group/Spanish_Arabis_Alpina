#!/bin/bash/

dir=$1
samples=$2
maxjob=$3
cd $dir

gatk2='PATH_TO_GATK/gatk-4.2.0.0/gatk'



mergebamalignment () {
        local j="$1"
        gunzip -k $j
        uj=$(echo $j | rev | cut -d '.' -f2- | rev)
	sample=$(echo $uj | rev | cut -d '_' -f2- | rev)
	$gatk2 --java-options "-XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:ActiveProcessorCount=1" MergeBamAlignment \
        	-R Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta \
        	--UNMAPPED_BAM ${uj} \
        	--ALIGNED_BAM ${sample}_samtofastq_bwamem.sam \
        	-O ${sample}_mergebamalignment.bam \
        	--CREATE_INDEX true \
        	--ADD_MATE_CIGAR true \
        	--CLIP_ADAPTERS false \
        	--CLIP_OVERLAPPING_READS true \
        	--INCLUDE_SECONDARY_ALIGNMENTS true \
        	--MAX_INSERTIONS_OR_DELETIONS -1 \
        	--PRIMARY_ALIGNMENT_STRATEGY BestMapq \
        	--ATTRIBUTES_TO_RETAIN XS \
		--TMP_DIR /scratch
	rm ${uj}
	gzip ${sample}_samtofastq_bwamem.sam
}



for i in $(eval ls ${samples}_fastqtosam.sam.gz)
do
        mergebamalignment $i &
        background=( $(jobs -p) )
        if (( ${#background[@]} == $maxjob )); then
                wait -n
        fi
done
wait

