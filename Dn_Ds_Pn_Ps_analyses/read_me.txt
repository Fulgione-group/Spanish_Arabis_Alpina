mont_a.py: this script converted the montbretiana fasta file to a dataframe
mont_a_to_vcf3.py: this script was used to add the montbretiana into the vcf.

mont_a_vcf.sh: this is the script to run snpeff.

ka_ks_7.py and ka_ks_7_gzip.py: These two scripts takes in the vcf and calculates the derived allele freq per site
this also tells if the site is non synonymous or synonymous as well as polymorphism or substitution.
one of these works on the gziped file other doesnot rest everything is same


ka_ks_7_cantabria_mont_a.sh and ka_ks_7_france_mont_a.sh: these are the bash scripts of how ka_ks_7.py was run on the cantabrian and french subset

dn_ds_2.py: this script uses the allele frequency file to calculate the Dn DS and Pn Ps per gene

dn_ds_2_france_mont_a_kakas7.sh, dn_ds_2_cantabria_mont_a_kakas7.sh: these are the bash scripts of how dn_ds_2.py was run on the cantabrian and french subset

cantabrians_vs_france_mont_a.ipynb: this is a jupyter notebook of the calculation of the dN/dS ratios

go_enrichment_pn_ps_paramask_filtered_top_go_final.R: this is script for GO enrichment analyses