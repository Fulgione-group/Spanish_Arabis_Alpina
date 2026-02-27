import pandas as pd
import numpy as np
import scipy
from scipy.stats import shapiro
from scipy.stats import norm
import re
from  scipy.stats import yeojohnson


output=open("pin_pis_per_gene_filter_dist_pvalue_stats_paramask_bootstrap_z.tsv","w")
output.write("gene"+"\t"+"mean"+"\t"+"std_dev"+"\t"+"z_p_value_gw"+"\t"+"z"+"\n")
output.close()
dist=open("pin_pis_dist_diff_france_spain_filter_bootstrap.txt","r")
dist=dist.read()
dist=dist.split("[")
dist=dist[1:]
#print(dist[0])
genes_df=pd.read_csv("pin_pis_stats_diff_france_spain_filter_bootstrap.tsv", sep="\t")

paramask=pd.read_csv("/netscratch/dep_coupland/grp_fulgione/mehak/data/spanish_project/pin_pis_bootstrap/paramask_filtered_genes.txt",sep="\t")
paramask.columns=["index","gene"]

genes_df=genes_df[genes_df["gene"].isin(list(paramask["gene"]))]
gw=np.nanmean(np.isfinite(genes_df["stat"]))
print(gw)
for i in genes_df.index:
    dist_gene=np.array(dist[i].rstrip().rstrip("]").split())
    dist_gene= np.array([float(x) for x in dist_gene])
    mu, sigma = norm.fit(dist_gene[np.isfinite(dist_gene)])
    z=(gw-mu)/sigma
    p_two = 2 * (1 - norm.cdf(abs(z)))
    output=open("pin_pis_per_gene_filter_dist_pvalue_stats_paramask_bootstrap_z.tsv","a")
    output.write(genes_df["gene"][i]+"\t"+str(mu)+"\t"+str(sigma)+"\t"+str(p_two)+"\t"+str(z)+"\n")
    output.close()
