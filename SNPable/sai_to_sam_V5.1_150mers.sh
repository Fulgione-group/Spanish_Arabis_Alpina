#!/bin/bash
bwa samse /netscratch/dep_coupland/grp_fulgione/bastiaan/data/reference/Alpina_V5.1/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta AaSplit150.sai AaSplit150.txt | bgzip -c > AaSplit150.sam.gz
