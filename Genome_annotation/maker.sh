#!/bin/bash
#BSUB -J maker
#BSUB -M 50000
#BSUB -W 72:00
#BSUB -n 12
#BSUB -R "rusage[mem=20000]"
#BSUB -q bigmem
#BSUB -o maker_%J.stdout
#BSUB -e maker_%J.stderr
#BSUB -N


ml maker
ml

ROUND=3
GENOME=pb
BASE=${GENOME}.ngenome
LOG=${BASE}_master_datastore_index.log


#export OMP_NUM_THREADS=2 
export LIBDIR=maker_v3.01.03/share/RepeatMasker/Libraries/
export REPEATMASKER_LIB_DIR=maker_v3.01.03/share/RepeatMasker/Libraries/
export REPEATMASKER_MATRICES_DIR=maker_v3.01.03/share/RepeatMasker/Matrices/
export PERL5LIB=perl5/lib/perl5/x86_64-linux-thread-multi:$PERL5LIB
export AUGUSTUS_CONFIG_PATH=augustus-config


maker maker_opts.ctl maker_bopts.ctl maker_exe.ctl -cpus 12 
