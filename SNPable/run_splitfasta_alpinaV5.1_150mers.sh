#!/bin/bash
cd /netscratch/dep_coupland/grp_fulgione/bastiaan/data/reference/SNPable
bsub -J "job_splitfasta" -o splifasta_out.txt -e S_splitfasta_err.txt -q normal -n 4 -R "rusage[mem=10240]" "/netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/seqbility-20091110/splitfa \
/netscratch/dep_coupland/grp_fulgione/bastiaan/data/reference/Alpina_V5.1/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta 150 > AaSplit150.txt"

