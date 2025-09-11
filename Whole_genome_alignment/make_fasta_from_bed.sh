#!/bin/bash

#grep chr length
grep 'chr' Arabis_alpina.MPIPZ.version_5.1.chr.all.genome.sizes


#
awk 'BEGIN{print ">Montbretiana alpina coordinate chr8"; FS="\t";p=0; cl=49795044}{diff=($2-p);if(diff!=0){for(i=1;i<=diff;i++){printf"%s","N"}}; \
split($4, a, "");for(i=1;i<=length(a);i++){printf"%s", a[i]};p=$3}END{diff=(cl-p);if(diff!=0){for(i=1; i<=diff;i++){printf"%s", "N"}}}' \
all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_${chr}_sorted.bed \
 > all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_${chr}_sorted.fasta


#check length

awk 'BEGIN{FS=""}{if(NR==2){print NF}}' all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_${chr}_sorted.fasta 
#or check with wc
tail +2 all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_${chr}_sorted.fasta | wc


#put 60 chars per line
awk 'BEGIN{FS=""}{if(NR==1){print $0}else{printf"%s", $1;for(i=2; i<=NF; i++){if(((i-1) % 60)==0){printf"\n"};printf"%s", $i}}}'\
 all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_${chr}_sorted.fasta \
> all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_${chr}_sorted_mlines.fasta


#check if character number still adds up
awk 'BEGIN{c=0; FS=""}{if(NR>1){c=(c+NF)}}END{print c}' all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_${chr}_sorted_mlines.fasta

#substitute gaps with Ns
sed 's/-/N/g' all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_${chr}_sorted_mlines.fasta  > all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_${chr}_sorted_mlines_noGaps.fasta


#check how many Ns are in the genome in total 25488358 missing sites out of 311642060 total sites in the chromosomes 1-8
awk 'BEGIN{FS="";c=0}{if(NR>1){for(i=1;i<=NF;i++){if($i=="N"){c++}}}}END{print c}' all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_chrall_sorted_mlines_noGaps.fasta

## count how many positions are in the reference
grep 'chr' Arabis_alpina.MPIPZ.version_5.1.chr.all.genome.sizes | awk 'BEGIN{c=0}{c=($2+c)}END{print c}'
