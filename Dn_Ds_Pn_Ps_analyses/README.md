# Scripts for dN/dS, Pn/Ps, and GO enrichment analyses

## Pre-processing
- **mont_a.py**  
  Converts the *A. montbretiana* FASTA file to a dataframe.  

- **mont_a_to_vcf3.py**  
  Adds the *A. montbretiana* sequence into the VCF.  

- **mont_a_vcf.sh**  
  Runs **SnpEff (v.5.2.a)** annotation on the VCF.

## Derived allele frequency and site classification
- **ka_ks_7.py** and **ka_ks_7_gzip.py**  
  Take the VCF as input and calculate the **derived allele frequency per site**.  
  They also classify each site as:  
  - **Non-synonymous or synonymous**  
  - **Polymorphism (Pn/Ps) or substitution (Dn/Ds)**  
  - `ka_ks_7.py`: works on uncompressed VCF.  
  - `ka_ks_7_gzip.py`: works on gzipped VCF.  

- **ka_ks_7_cantabria_mont_a.sh** and **ka_ks_7_france_mont_a.sh**  
  Bash scripts to run `ka_ks_7.py` on the Cantabrian and French subsets, respectively.

## Gene-level dN/dS and Pn/Ps calculations
- **dn_ds_2.py**  
  Uses the allele frequency file to calculate **Dn, Ds, Pn, and Ps per gene**.  

- **dn_ds_2_france_mont_a_kakas7.sh** and **dn_ds_2_cantabria_mont_a_kakas7.sh**  
  Bash scripts to run `dn_ds_2.py` on the Cantabrian and French subsets.

- **cantabrians_vs_france_mont_a.ipynb**  
  Jupyter notebook to calculate and visualize **dN/dS** and **Pn/Ps** ratios between Cantabrian and French populations.

## GO enrichment
- **go_enrichment_pn_ps_paramask_filtered_top_go_final.R**  
  Performs **GO enrichment analysis** based on the **Pn/Ps results**.
