# Spanish_Arabis_Alpina

Scripts for all analyses supporting the manuscript  
*“Adaptation to seasonal drought in Arabis alpina is linked to the demographic history and climatic changes since the last glacial maximum”*

Each bottom-level subdirectory includes a README file.

Scripts are included for:

1. **SNP_calling**
   - SNP calling with GATK

2. **Phasing**
   - Two-step phasing approach using WhatsHap (read-pair information) and SHAPEIT (HMM-based phasing and imputation)

3. **SNPable (mappability mask)**
   - K-mer unique mapping based approach

4. **ParaMask**
   - Detection of multicopy regions

5. **3PCLR_Analysis**
   - Detection of signatures of positive selection

6. **Demography**
   - Based on genome-wide genealogies using Relate

7. **Structure_Analysis**
   - Diversity estimates, Fis estimates, PCoA, and Neighbor-Joining tree based on 1-IBS distances

8. **Dn_Ds_Pn_Ps_analyses**
   - Differential synonymous vs. non-synonymous polymorphism and divergence analysis

9. **Genotype_association_analysis (NAC055 and FRL1 genes)**
   - Linear modeling using structure covariates to test associations of genes with monthly aridity and flowering time

