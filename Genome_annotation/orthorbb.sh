#!/bin/bash
#BSUB -J orthorbb
#BSUB -M 50000
#BSUB -W 24:00
#BSUB -n 16
#BSUB -R "rusage[mem=20000]"
#BSUB -q bigmem
#BSUB -o orthorbb_%J.stdout
#BSUB -e orthorbb_%J.stderr
#BSUB -N


cd Arabis_alpina.MPIPZ.version_5.1.chr.all.maker.output
pwd -P 

INPUT=Aa_5.1_helixer.maker.all.maker.proteins.filtered.merged.mod.fasta

UNIPROT=uniprot-reviewed_3A-2022.08.10-15.43.36.82.fasta

orthorbb -p $INPUT -a $UNIPROT -q Aa_proteins -r UniProt -t 16

sed -s 's/|[OR]BB\s/\t/g' Aa_proteins_UniProt_homology_annotation.fmt6.txt > Aa_proteins_UniProt_homology_annotation.fmt6.mod.txt

OUTPUT=Aa_proteins_UniProt_homology_annotation.fmt6.mod.txt
GFF=Aa_5.1_helixer.maker.genes.tidy.v2.curated.sorted.mod.gff.ipr
PROT=Aa_5.1_helixer.maker.all.maker.proteins.filtered.merged.mod.fasta
TRANS=Aa_5.1_helixer.maker.all.maker.transcripts.filtered.merged.mod.fasta

maker_functional_gff $UNIPROT $OUTPUT $GFF > $GFF.ortho
maker_functional_fasta $UNIPROT $OUTPUT $PROT > $PROT.ortho
maker_functional_fasta $UNIPROT $OUTPUT $TRANS > $TRANS.ortho
