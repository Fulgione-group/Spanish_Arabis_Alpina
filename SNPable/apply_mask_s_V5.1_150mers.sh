#!/bin/bash
for i in $(ls mask_150_r0.*)
do
	/netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/seqbility-20091110/apply_mask_s $i /netscratch/dep_coupland/grp_fulgione/bastiaan/data/reference/Alpina_V5.1/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta > $(echo $i | cut -d'_' -f3)_150_genome_mask.fa
done
