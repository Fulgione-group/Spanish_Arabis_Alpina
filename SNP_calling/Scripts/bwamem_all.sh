#!/bin/bash/

dir=$1
echo "dir"
samples=$2
maxjob=$(expr $(expr $3 - 1 ) / 8 )
cd $dir

gatk2='PATH_TO_GATK/gatk-4.2.0.0/gatk'


bwamap () {
        local j="$1"
        gunzip -k $j
        uj=$(echo ${j} | rev | cut -d '.' -f2- | rev)
        sample=$(echo ${uj} | rev | cut -d '_' -f2- | rev)
        bwa index ${uj} &&
        bwa mem -M -t 8 -p Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta \
        ${uj} > ${sample}_bwamem.sam 2> ${sample}_bwamem_stderr.txt
	rm ${uj}
	rm ${uj}.ann
	rm ${uj}.amb
	rm ${uj}.bwt
	rm ${uj}.pac
	rm ${uj}.sa
}


for i in $(eval ls ${samples}_samtofastq_interleaved.fastq.gz)
do
        bwamap $i &
        background=( $(jobs -p) )
        if (( ${#background[@]} == $maxjob )); then
                        wait -n
        fi
done
wait



