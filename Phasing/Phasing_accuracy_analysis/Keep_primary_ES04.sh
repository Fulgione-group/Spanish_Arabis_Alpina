#!/bin/bash
conda activate minimap2_env

samtools view -b -F 0x900 /phasing_accuracy/mapped_reads/ERR15664662_ES04.bam\
  | samtools sort -o /phasing_accuracy/mapped_reads/ERR15664662_primary_sorted_ES04.bam

samtools index /phasing_accuracy/mapped_reads/ERR15664662_primary_sorted_ES04.bam
