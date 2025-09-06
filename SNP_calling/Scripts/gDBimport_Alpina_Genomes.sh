#!/bin/bash

export TILEDB_DISABLE_FILE_LOCKING=1

cd /netscratch/dep_coupland/grp_fulgione/bastiaan/data/GATK_DBs/
gatk2='PATH_TO_GATK/gatk-4.2.0.0/gatk'

$gatk2 --java-options "-XX:ConcGCThreads=1 -XX:ActiveProcessorCount=3 -Xmx150G" GenomicsDBImport \
	--sample-name-map gDB_sample_map_1000Genomes_sorted.txt \
	--max-num-intervals-to-import-in-parallel 1 \
	-L chr1 \
	-L chr2 \
	-L chr3 \
	-L chr4 \
	-L chr5 \
	-L chr6 \
	-L chr7 \
	-L chr8 \
	--genomicsdb-workspace-path PATH_TO_GATK_DBs/AlpinaGenomes \
	--reference Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta


