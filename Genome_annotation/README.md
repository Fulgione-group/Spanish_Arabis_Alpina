# Collection of scripts for the genome annotation

Structural and functional genome annotation of *Arabis alpina* following these steps:


## Repeat masking
- This aids the identification of genes by reducing the search space in the genome and provides a valuable resource on its own. A de novo repeat identification using `NCBI BLAST` and the `BuildDatabase` command was run. These repeat sequences are used to mask the assembly.
`repeatmasker.sh`

## Structural annotation
`MAKER` is an extensive pipeline that bundles a series of tools to perform genome annotation. It includes evidence-based and ab initio driven approaches to predict genes and other genomic features. Its output is the structural annotation of a genome in standard file formats, which can be further analysed to add functional annotations such as the process a gene is involved in, its molecular function, or location of expression. 

- Repeat masking was performed outside of the Maker pipeline and the resulting GFF was supplied (see above).

- Prediction of gene structure using `Helixer`, based on Deep Neural Networks:
`helixer.sh`

- Ab initio gene prediction using `Augustus`, based on Hidden-Markov Models: 
`augustus.sh`

- Control files needed to run `Maker`:
`maker_bopts.ctl`
`maker_evm.ctl`
`maker_opts.ctl`

- Script to start the actual `Maker` pipeline:
`maker.sh`

## Functional annotation

- The Maker output was further processed with accessory scripts from Maker as well as various functions from the `AGAT` suite.
postprocessing.sh

- Classification of protein families:
interproscan.sh

- Inferring protein homology:
orthorbb.sh
