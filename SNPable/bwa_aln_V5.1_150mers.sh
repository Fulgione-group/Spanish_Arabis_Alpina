#!/bin/bash
cd /netscratch/dep_coupland/grp_fulgione/bastiaan/data/reference/SNPable
bwa aln -t 20 -R 1000000 -O 3 -E 3 /netscratch/dep_coupland/grp_fulgione/bastiaan/data/reference/Alpina_V5.1/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta AaSplit150.txt > AaSplit150.sai
