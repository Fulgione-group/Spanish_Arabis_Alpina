#!/bin/bash

dir=$1
echo "$dir"
cd $dir
structure_file=$2
maxjob=$(expr $3 - 1 )


gatk2='PATH_TO_GATK/gatk-4.2.0.0/gatk'



markduplicated () {
        eval j="$1"
	eval Out="$2"
        $gatk2 --java-options "-XX:ConcGCThreads=2 -XX:ActiveProcessorCount=4" MarkDuplicates \
        	$j \
        	-O ${Out}_dm.bam \
        	-M ${Out}_dm_metrics.txt \
        	--TMP_DIR /scratch
}


counter=1
filecounter=0
declare -a filearray=();
while read -r line;
do
        if [ "$counter" -eq 1 ] || [[ "$SM" == "$(echo $line| awk '{print $1}')" ]]
        then
                filearray[$filecounter]="-I $(echo $line| awk '{print $2}')"
                filecounter=$((filecounter + 1))
                SM=$(echo $line| awk '{print $1}')
        else
        	markduplicated "\${filearray[@]}" "\${SM}" &
        	background=( $(jobs -p) )
        	if (( ${#background[@]} == $maxjob )); then
                	        wait -n
        	fi
		declare -a filearray=()
                filecounter=0
                SM=$(echo $line| awk '{print $1}')
                filearray[$filecounter]="-I $(echo $line| awk '{print $2}')"
                filecounter=$((filecounter + 1))
        fi
        counter=$((counter+1))
done < $structure_file
wait

#process last sample outside loop
markduplicated "\${filearray[@]}" "\${SM}"
background=( $(jobs -p) )
if (( ${#background[@]} == $maxjob )); then
	wait -n
fi
wait
