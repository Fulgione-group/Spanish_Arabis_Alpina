#!/bin/bash/

dir=$1
outdir=$2
samples=$3
maxjob=$4
Platform=$5
Company=$6
OutID=$7
cd $dir

gatk2='PATH_TO_GATK/gatk-4.2.0.0/gatk'

fastqtosam_illumina () {
		local j="$1"
		fastq=$(echo $j | rev | cut -d '.' -f2- | rev )
		if [[ $fastq =~ (.+)R1(.*) ]]
        	then
                	fastq2=$(echo $fastq | sed 's/R1/R2/g')
        	else
                	r1=FALSE
        	fi
        	if [[ $fastq =~ (.+)fwd(.*) ]]
        	then
                	fastq2=$(echo $fastq | sed 's/fwd/rev/g')
        	else
                	fwd=FALSE
        	fi
        	if [ r1 == FALSE ] && [ fwd == FALSE]
        	then
                	echo "error: readgroup string not found"
                	break
        	fi
		SN=$(echo $fastq | cut -d '_' -f1,2)
		#if [ -z ${OutID+x} ] # didnt work because empty string was detected with +x
		if [ -z ${OutID} ]
		then
			OUT=${outdir}/$(echo $fastq | cut -d '_' -f1,2)_fastqtosam.sam
		else
			OUT=${outdir}/$(echo $fastq | cut -d '_' -f1,2)_${OutID}_fastqtosam.sam
		fi
		RG=$(head -1 $fastq.gz | awk '{split($0,a,"\:"); split(a[3],b,""); printf "%s" ,b[1];for(i=2;i<=length(b);i++){if(b[i]~/[A-Z]/){break}else{printf "%s",b[i]}};printf ".%s\n" ,a[4]}' 2> /dev/null)
		echo "f1=$fastq" > ${OUT}_log_file.txt
		echo "f2=$fastq2" >> ${OUT}_log_file.txt
		echo "RG=$RG" >> ${OUT}_log_file.txt
		echo "Out=$OUT" >> ${OUT}_log_file.txt
		echo "SN=$SN" >> ${OUT}_log_file.txt
		echo "Company=$Company" >> ${OUT}_log_file.txt
		$gatk2 --java-options "-XX:ConcGCThreads=1 -XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -Xmx6G" FastqToSam \
			-F1 $fastq.gz \
	   		-F2 $fastq2.gz \
	   		--OUTPUT $OUT \
	   		--READ_GROUP_NAME $RG \
	   		--SAMPLE_NAME $SN \
 	  		--LIBRARY_NAME Truseq-$SN \
 	   		--PLATFORM Illumina \
 	   		--SEQUENCING_CENTER $Company
		gzip $OUT
}

fastqtosam_BGI () {
		local j="$1"
   	    	fastq=$(ls ${j}/*_1.fq.gz | sed 's/.gz//g')
        	fastq2=$(ls ${j}/*_2.fq.gz | sed 's/.gz//g')
        	SN=$j
        	if [ -z ${OutID+x} ]
        	then
                	OUT=${outdir}/${SN}_fastqtosam.sam
        	else
                	OUT=${outdir}/${SN}_${OutID}_fastqtosam.sam
        	fi
       		RG=$(gzip -dc $fastq.gz | head -1 | awk '{split($0,a, "C[[:digit:]]{3}"); print a[1]}')
        	echo "f1=$fastq" > ${OUT}_log_file.txt
        	echo "f2=$fastq2" >> ${OUT}_log_file.txt
        	echo "RG=$RG" >> ${OUT}_log_file.txt
        	echo "Out=$OUT" >> ${OUT}_log_file.txt
		echo "SN=$SN" >> ${OUT}_log_file.txt
		echo "Company=$Company" >> ${OUT}_log_file.txt
        	$gatk2 --java-options "-XX:ConcGCThreads=1 -XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -Xmx6G" FastqToSam \
           		-F1 $fastq.gz \
           		-F2 $fastq2.gz \
           		--OUTPUT $OUT \
           		--READ_GROUP_NAME $RG \
           		--SAMPLE_NAME $SN \
           		--LIBRARY_NAME MGISEQ200-$SN \
		        --SEQUENCING_CENTER $Company
                gzip $OUT
}
if [[ $Platform == "Illumina" ]]
then
	for i in $(ls ${samples}*.f*q.gz)
	do
		fastqtosam_illumina "$i" &
               	background=( $(jobs -p) )
                if (( ${#background[@]} == $maxjob )); then
                        wait -n
                fi
	done
	wait
elif [[ $Platform == "DNBseq" ]]
then
	for i in $( ls -d ${samples})
	do
		fastqtosam_BGI "$i" &
	      	background=( $(jobs -p) )
	        if (( ${#background[@]} == $maxjob )); then
        	        wait -n
        	fi
	done
	wait
else
	echo "error no Platform detected"
fi
