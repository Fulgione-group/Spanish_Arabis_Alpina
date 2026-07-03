# XP-EHH analysis

This folder contains scripts to process phased SNPs for XP-EHH selection scans using **selscan** (Szpiech & Hernandez, 2014).

## Prepare input

Randomly select one haplotype per individual to avoid biases caused by recent inbreeding during plant propagation:

- `haploidize.sh`

Recode the haploidized genotypes into diploid format for **selscan**:

- `make_diploid_coded.sh`

Generate map files using recombination maps:

- `generate_maps.sh`

Split VCFs by chromosome:

- `Split_by_chr.sh`

## Run XP-EHH

Run XP-EHH using two VCF files and one recombination map. The option `--max-gap 500000` prevents EHH integration across long physical gaps, while `--trunc-ok` allows truncated integration at chromosome ends.

- `run_xpehh_parallel_GROUP1_vs_GROUP2.sh`

## Analyze overlap with 3P-CLR peaks

Load the XP-EHH results and identify overlap with the highest-scoring 3P-CLR peaks:

- `XPEHH_overlap_with_high_3PCLR_peaks.R`
