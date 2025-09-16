#!/bin/bash
#BSUB -J RepeatMasker
#BSUB -M 50000
#BSUB -W 72:00
#BSUB -n 4
#BSUB -R "rusage[mem=20000]"
#BSUB -q bigmem
#BSUB -o RepeatMasker_%J.stdout
#BSUB -e RepeatMasker_%J.stderr
#BSUB -N

ml maker
ml

export LIBDIR=maker_v3.01.03/share/RepeatMasker/Libraries/
export REPEATMASKER_LIB_DIR=maker_v3.01.03/share/RepeatMasker/Libraries/
export REPEATMASKER_MATRICES_DIR=maker_v3.01.03/share/RepeatMasker/Matrices/
export PERL5LIB=perl5/lib/perl5/x86_64-linux-thread-multi/:$PERL5LIB

FASTA=Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta


LIB=Aalpina-families.fa
DIR=

BuildDatabase -name Aalpina -engine rmblast $FASTA
RepeatModeler -database Aalpina -engine ncbi -pa 4

RepeatMasker -e ncbi -gff -pa 4 -lib $LIB -s -dir $DIR -xsmall $FASTA

# script downloaded from:
# https://github.com/darencard/GenomeAnnotation/blob/8b0e5cb5b887e1d576e8483ded8e5cc754e75d23/rmOutToGFF3custom

rmOutToGFF3.pl Arabis_alpina.test.10MB.fasta.out > Arabis_alpina.test.10MB.fasta.gff3
cat Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gff3 | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.reformat.gff3

