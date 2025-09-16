#!/bin/bash
#BSUB -J augustus
#BSUB -M 50000
#BSUB -W 24:00
#BSUB -n 8
#BSUB -R "rusage[mem=20000]"
#BSUB -q bigmem
#BSUB -o augustus_%J.stdout
#BSUB -e augustus_%J.stderr
#BSUB -N

ml busco/v5.1.2
ml

BASE=Aalpina
DIR=augustus
GFF=${BASE}.maker.noseq.gff 

mkdir -p $DIR
cd $DIR
pwd -P

FASTA=${BASE}.fasta 

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' $GFF | awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | bedtools sort -chrThenSizeA -i - | uniq | bedtools getfasta -fi $FASTA -bed - -fo ${BASE}.maker.transcripts1000.fasta

export AUGUSTUS_CONFIG_PATH=augustus-config

busco -i ${BASE}.maker.transcripts1000.fasta -l brassicales_odb10 -o ${BASE}_maker -m genome --augustus --augustus_parameters='--progress=true' --augustus_species arabidopsis -c 8 --long --force

cd ${BASE}_maker/run_brassicales_odb10/augustus_output/retraining_parameters/BUSCO_${BASE}_maker/

echo "renaming BUSCO output files"
rename "s/BUSCO_${BASE}_maker/${BASE}/g" *
sed -i "s/BUSCO_${BASE}_maker/${BASE}/g" ${BASE}_parameters.cfg
sed -i "s/BUSCO_${BASE}_maker/${BASE}/g" ${BASE}_parameters.cfg.orig1

NEWDIR=
mkdir -p $NEWDIR

cp ${BASE}_${RUN}_rnd${ROUND}* $NEWDIR