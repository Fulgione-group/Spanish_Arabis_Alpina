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

## Logistic Regression(GEE)
- **dn_ds_logistic_simple_input_format.py** and **pn_ps_logistic_simple_input_format.py**  
  Filters for SNPs and converts the data to logistic regression input format for substitutions and polymorphisms.  

- **dn_ds_logistic_simple_2_g2_g5_filter_for_low_pvt_var_5_pr_wald_0.016631363015617805_gee.py** and **pn_ps_logistic_simple_2_g2_g10_filter_for_low_pvt_var_10_wald0.4_gee.py**  
  Scripts to perform GEE logistic regression on substitutions and polymorphisms

- **log_reg_dn_ds_gee.ipynb** and **log_reg_pn_ps_gee.ipynb**  
  Jupyter notebook to visualize the coefficient distribution and do multiple hypothesis correction.


