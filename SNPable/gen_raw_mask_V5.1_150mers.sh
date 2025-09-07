#!/bin/bash

gzip -dc AaSplit150.sam.gz | /netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/seqbility-20091110/gen_raw_mask.pl > rawMask_150_V5.1.fa
