import pandas as pd
import numpy as np
from scipy.stats import bootstrap

output_dist=open("pin_pis_dist_diff_france_spain_filter_bootstrap_filtered_v2.txt","w")
output=open("pin_pis_stats_diff_france_spain_filter_bootstrap_filtered_v2.tsv","w")
output.write("gene"+"\t"+"stat"+"\t"+"std_err"+"\t"+"ci"+"\n")
output_dist.close()
def my_stat(pin_france, pis_france,pin_cant,pis_cant):
    pin_pis_france=np.mean(pin_france)/np.mean(pis_france)
    pin_pis_cant=np.mean(pin_cant)/np.mean(pis_cant)
    diff=pin_pis_france-pin_pis_cant
    return diff


df_france=pd.read_csv("pi_per_site3_france_mont_a.txt", sep="\t", header=None )
df_cantabria=pd.read_csv("pi_per_site3_cantabria_mont_a.txt", sep="\t", header=None )
df_france.columns=["chr","pos","pi_france","gene","ns"]
df_cantabria.columns=["chr","pos","pi_cant","gene","ns"]
paramask=pd.read_csv("paramask_filtered_genes.txt",sep="\t")
paramask.columns=["index","gene"]
genes=list(paramask["gene"])
df_france["pi_france"] = pd.to_numeric(df_france["pi_france"],errors="coerce")
df_cantabria["pi_cant"] = pd.to_numeric(df_cantabria["pi_cant"],errors="coerce")
df_france=df_france.dropna()
df_cantabria=df_cantabria.dropna()
for gene in genes:
    df_france_gene=df_france[df_france["gene"]==gene]
    df_cantabria_gene=df_cantabria[df_cantabria["gene"]==gene]
    if len(df_france_gene[df_france_gene["pi_france"]!=0])>2 and len(df_cantabria_gene[df_cantabria_gene["pi_cant"]!=0])>2:
        df_f_c_gene=pd.merge(df_france_gene, df_cantabria_gene, on="pos")
        df_f_c_gene_n=df_f_c_gene[df_f_c_gene["ns_x"]=="n"]
        df_f_c_gene_s=df_f_c_gene[df_f_c_gene["ns_x"]=="s"]
        pin_france_np=df_f_c_gene_n["pi_france"].dropna().to_numpy()
        pis_france_np=df_f_c_gene_s["pi_france"].dropna().to_numpy()
        pin_cant_np=df_f_c_gene_n["pi_cant"].dropna().to_numpy()
        pis_cant_np=df_f_c_gene_s["pi_cant"].dropna().to_numpy()
        print(pin_france_np.shape, pis_france_np.shape, pin_cant_np.shape, pis_cant_np.shape)
        if pin_france_np.shape[0]>=5 and pis_france_np.shape[0]>=5 and pin_cant_np.shape[0]>=5 and pis_cant_np.shape[0]>=5:
            result=bootstrap((pin_france_np, pis_france_np,pin_cant_np,pis_cant_np,), statistic=my_stat,n_resamples= 100)
            output.write(gene+"\t"+str(my_stat(pin_france_np, pis_france_np,pin_cant_np,pis_cant_np))+"\t"+str(result.standard_error)+"\t"+str(result.confidence_interval)+"\n")
            output_dist=open("pin_pis_dist_diff_france_spain_filter_bootstrap_filtered_v2.txt","a")
            output_dist.write(str(result.bootstrap_distribution)+"\n")
            output_dist.close()
output.close()
