#!/bin/bash
from="$1"
to="$2"

run_whatshap () {
	local j="$1"
	#read bams from a text file
	readarray -t arguments < $PATH_TO_WD/BamList_sample${j}.txt
	$PATH_TO_WHATSHAP/whatshap phase\
		-o $PATH_TO_WD/GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.p${j}.prephased.b.vcf.gz\
 		--reference $PATH_TO_REFERENCE/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta\
 		--tag=PS\
		$PATH_TO_WD/GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.p${j}.b.vcf.gz\
		"${arguments[@]}"
}

for (( i=$from; i<=$to; i++ ))
do
        run_whatshap $i &
        background=( $(jobs -p) )
        if (( ${#background[@]} == 10 )); then
                wait -n
        fi
done
wait

