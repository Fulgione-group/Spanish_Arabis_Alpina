# Overview 
We used GATK [v4.2.0.0](https://gatk.broadinstitute.org/hc/en-us/sections/360012354372-4-2-0-0) ([DePristo et al., 2011](https://www.nature.com/articles/ng.806)) for SNP calling in the two *Arabis alpina* populations. Read groups were assigned, adapters were soft clipped, duplicates were marked and reads were aligned to the pajares reference version [v5.1](http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz) ([W.-B. Jiao and K. Schneeberger, 2017](https://www.sciencedirect.com/science/article/pii/S1369526616301315?via%3Dihub)) ) using the mem algorithm of bwa v0.7.17 ([H. Li and R. Durbin, 2009](https://pubmed.ncbi.nlm.nih.gov/19451168/)). Additionally, we hard filter 

<br>
<br>

# Pipeline:
## fastqtosam_all.sh

Script convert fastq files to unmapped sam files and assign read groups from header information (valid for DNBseq headers 2020). 
<br>
## markadapters_all.sh and discount_all.sh
\*Note: Adapters were already cleaned by the sequencing provider, Beijing Genomics Institute (BGI). 
Scripts to markadapters and soft clip them, supplied with five prime and three prime adapters from DNBseq.
<br>
## bwamem_all.sh
Script to map samples to the *A. alpina* pajares reference.
<br>
## mergebamalignment_all.sh and sortmabam_all.sh
Script to merge mapped bam file and raw unmapped sam files and sort the merged bam file
<br>
## markadapters
Mark PCR duplicates (reads from the same DNA sequencing freqment)
<br>
## haplocallergvcf_all2.sh
Per sample genotyping using GATK's haplotype caller with "--output-mode EMIT_ALL_ACTIVE_SITES" to emit an all sites VCF (including invariant sites). Outputs in raw genotype calls in GVCFs
<br>
## genotypeGVCF_Alpina_Genomes.sh
Creates a GATK genomic Database
<br>
## genotypeGVCF_Alpina_Genomes.sh
Performs joint genotyping for all samples in gDB created in the previous step
<br>
## run_FilterVCF_Alpina_Genomes.sh
Filter with custom scripts per genotype and per site: Genotype quality >=30, depth >=5, min-genome-wide-cov >= 10x, biallic SNps only
