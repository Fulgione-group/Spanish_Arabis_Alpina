#!/bin/bash

dir=$1
echo "dir"
cd $dir
samples=$2
gatk2='PATH_TO_GATK/gatk-4.2.1.0/gatk'


$gatk2 --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -XX:ConcGCThreads=2 -XX:ActiveProcessorCount=4 -Xms10G -Xmx80G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
    -R Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta \
    -V gendb://PATH_TO_GATK_DBs/AlpinaGenomes \
    --use-new-qual-calculator \
    -O GATK4.2_Alpina_Genomes.vcf \
    -all-sites\
    --tmp-dir /scratch\
    -L chr1\
    -L chr2\
    -L chr3\
    -L chr4\
    -L chr5\
    -L chr6\
    -L chr7\
    -L chr8
