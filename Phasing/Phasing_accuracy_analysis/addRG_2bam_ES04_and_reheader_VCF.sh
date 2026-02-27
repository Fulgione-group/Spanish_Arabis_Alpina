#!/bin/bash

conda activate env_bcftools

cd /phasing_accuracy/mapped_reads/


samtools addreplacerg \
  -r 'ID:ES04_hifi\tSM:ES04_hifi\tPL:PACBIO\tLB:hifi\tPU:unit1' \
  -o ERR15664662_primary_sorted_ES04.RG.bam \
  ERR15664662_primary_sorted_ES04.bam

samtools index ERR15664662_primary_sorted_ES04.RG.bam

