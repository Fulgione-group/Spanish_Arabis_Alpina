#!/bin/bash
#BSUB -J interproscan
#BSUB -M 50000
#BSUB -W 12:00
#BSUB -n 16
#BSUB -R "rusage[mem=20000]"
#BSUB -q bigmem
#BSUB -o interproscan_%J.stdout
#BSUB -e interproscan_%J.stderr
#BSUB -N


ml interproscan
ml

PROT=Aa_5.1_helixer.maker.all.maker.proteins.filtered.merged.mod.fasta

interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -cpu 16 -i $PROT -o ${PROT}.iprscan

GFF=Aa_5.1_helixer.maker.genes.tidy.v2.curated.sorted.mod.gff
PROT=Aa_5.1_helixer.maker.all.maker.proteins.filtered.merged.mod.fasta.iprscan

ipr_update_gff $GFF ${PROT} > ${GFF}.ipr
iprscan2gff3 ${PROT} ${GFF}.ipr > ${GFF}.iprscan

quality_filter.pl -d ${GFF}_ipr > ${GFF}_ipr_def
quality_filter.pl -s ${GFF}_ipr > ${GFF}_ipr_std