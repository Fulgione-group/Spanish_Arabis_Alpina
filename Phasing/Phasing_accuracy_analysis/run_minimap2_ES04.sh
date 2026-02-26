#!/bin/bash
conda activate minimap2_env

minimap2 -t 16 -ax map-hifi Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta \
 /phasing_accuracy/longreads_raw/ERR15664661.fastq.gz \
  | samtools sort -@ 8 -o /phasing_accuracy/mapped_reads/ERR15664662_ES04.bam
samtools index /phasing_accuracy/mapped_reads/ERR15664662_ES04.bam

