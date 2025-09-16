#!/bin/bash
#BSUB -J makerPost
#BSUB -M 10000
#BSUB -W 12:00
#BSUB -n 4
#BSUB -R "rusage[mem=5000]"
#BSUB -q bigmem
#BSUB -o makerPost_%J.stdout
#BSUB -e makerPost_%J.stderr
#BSUB -N


ml maker
ml

RUN=
ROUND=
BASE=pb.ngenome
LOG=${BASE}_master_datastore_index.log

DIR=
cd $DIR
pwd -P

echo "postprocessing of output"

cd ${BASE}.maker.output

echo "gff3 and fasta merging"

gff3_merge -s -n -d $LOG > ${BASE}_rnd${ROUND}.maker.noseq.gff 
gff3_merge -g -s -n -d $LOG > ${BASE}_rnd${ROUND}.maker.genes.gff
fasta_merge -d $LOG

echo "evidence specific gff"

GFF=${BASE}_rnd${ROUND}.maker.noseq.gff 
# transcript alignments
awk '{ if ($2 ~ "est2genome") print $0 }' $GFF > ${BASE}_rnd${ROUND}.maker.est2genome.gff
# protein alignments
awk '{ if ($2 ~ "protein2genome") print $0 }' $GFF > ${BASE}_rnd${ROUND}.maker.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' $GFF > ${BASE}_rnd${ROUND}.maker.repeats.gff

VAR=${BASE}_rnd${ROUND}.maker.genes.gff
echo "AED"
AED_cdf_generator.pl -b 0.025 $VAR

echo "genestats"
echo "total genes: "
cut -f 3 $VAR | grep -c "gene"

# script downloaded from:
# https://github.com/darencard/GenomeAnnotation/blob/8b0e5cb5b887e1d576e8483ded8e5cc754e75d23/genestats

genestats $GFF > gene.stats

cat gene.stats | \
awk -v OFS="\t" '{ exon_sum += $3; exon_len += $4; intron_sum += $5; intron_len += $6 ; gene_len += $2 } END { print "total genes: " NR, "mean gene lengths: " gene_len / NR, "total exons: " exon_sum, "mean ex length: " exon_len / exon_sum, "total introns: " intron_sum, "mean in lengths: " intron_len / intron_sum }'
