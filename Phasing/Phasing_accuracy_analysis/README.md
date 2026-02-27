# General strategy
We estimated phasing accuracy by comparing heterozygotes sites of phased short reads (Read based with Whatshap and HMM based with shapeit) to read based phasing of 2 long read individuals


## Pipeline for one example accession ES04

### Mapping
- Align long reads to the reference using minimap2, sort and index bam: 
`run_minimap2_ES04.sh`

- Extract primary alignments only, sort and index bam:
`Keep_primary_ES04.sh`

### SNP calling
- call SNPs with bcftools (GATK is not compatible with minimap2 alignment scores):
`call_with_bcftools_ES04.sh`

- Filter sites with less than 30 genotype quality or less depth 5. (Note this is not exactly equivalent to GATKs genotype quality): 
`filter_GQ30DP5_bcftools_ES04.sh`

- reheader the VCF to give the sample a meaningful name:
`reheader_VCF_ES04.sh`


### run read group based phasing for longreads
- add proper Read group to bam for whatshap later
`addRG_2bam_ES04_and_reheader_VCF.sh`

- Read based phasing of the long-read VCF:
`whatshap_ES04_hifi.sh`


### filtering  before comparison
- retain only heterozygotes:
`retain_hets_phased_bcftools_ES04.sh`

- filter out SNPapble ParaMask regions:
`filter_bedtools_ParaMask_nogenes0kb_CAN_FR_ES04.sh`

- polarize VCF to montbretiana: 
`polarize_hifi_ES04.sh`
`add_anc_flag_header.sh`

### Short read VCF handling 
- extract het sites of the accession from a phased VCF:
`extract_ES04_shortreads.sh`


### Find overlap between short- and longreads
- Find an extract positions with bcftools:
`intersect_long_short_reads_ES04.sh`


### Run Whatshapp COmpare:
- finally run whatshap and look at switch rates
`whatshap_compare_ES04.sh`
