import pandas as pd
df_france_all_mut=pd.read_csv("output_kaks_7_france_mont_a.txt", sep="\t", header=None)
df_france_all_mut.columns=["chrom","pos","non_syn_allele_freq","syn_allele_freq","gene","syn_non_syn","polymorphism_substitution"]
df_cantabria_all_mut=pd.read_csv("output_kaks_7_cantabria_mont_a.txt", sep="\t", header=None)
df_cantabria_all_mut.columns=["chrom","pos","non_syn_allele_freq","syn_allele_freq","gene","syn_non_syn","polymorphism_substitution"]
df_cantabria_all_mut_p=df_cantabria_all_mut[df_cantabria_all_mut["polymorphism_substitution"]=="s"]
df_france_all_mut_p=df_france_all_mut[df_france_all_mut["polymorphism_substitution"]=="s"]
## we dont need two columns for allele freq we can easlily use 1
df_france_all_mut_p["daf"]=df_france_all_mut_p["non_syn_allele_freq"]+df_france_all_mut_p["syn_allele_freq"]
df_cantabria_all_mut_p["daf"]=df_cantabria_all_mut_p["non_syn_allele_freq"]+df_cantabria_all_mut_p["syn_allele_freq"]
df_cantabria_all_mut_p=df_cantabria_all_mut_p[["chrom","pos","gene","syn_non_syn","daf"]]
df_france_all_mut_p=df_france_all_mut_p[["chrom","pos","gene","syn_non_syn","daf"]]
df_france_all_mut_p.columns=['chrom', 'pos', 'gene', 'syn_non_syn', 'daf_france']
df_cantabria_all_mut_p.columns=['chrom', 'pos', 'gene', 'syn_non_syn', 'daf_cantabria']
df_cantabria_all_mut_p[df_cantabria_all_mut_p["daf_cantabria"]!=0]
df_france_all_mut_p[df_france_all_mut_p["daf_france"]!=0]
df=pd.merge(df_france_all_mut_p,df_cantabria_all_mut_p, how="outer", on=["chrom","pos"])
del df_france_all_mut
del df_cantabria_all_mut
del df_france_all_mut_p
del df_cantabria_all_mut_p
##removing the polymorphisms with 0 daf in both france and cant
df=df[(df["daf_france"]!=0) | (df["daf_cantabria"]!=0)]
df=df.fillna(0)
genes=df["gene_x"].unique()
output=open("output_log_reg_df_substitution.tsv","w")
output.write("chrom" + "\t" + "pos" + "\t"+"gene"+"\t"+"set"+"\t"+"ns"+"\n")
for gene in genes:
    df_gene=df[df["gene_x"]==gene]
    df_gene_f=df_gene[df_gene["daf_france"]!=0]
    for i in df_gene_f.index:
        output.write(df_gene_f["chrom"][i] + "\t" + str(df_gene_f["pos"][i])+ "\t"+df_gene_f["gene_x"][i] + "\t" +"france"+"\t"+df_gene_f["syn_non_syn_x"][i]+"\n")
        print(df_gene_f["chrom"][i] + "\t" + str(df_gene_f["pos"][i])+ "\t"+df_gene_f["gene_x"][i] + "\t" +"france"+"\t"+df_gene_f["syn_non_syn_x"][i])
    df_gene_c=df_gene[df_gene["daf_cantabria"]!=0]
    for i in df_gene_c.index:
        output.write(df_gene_c["chrom"][i] + "\t" + str(df_gene_c["pos"][i])+"\t"+ df_gene_c["gene_y"][i] + "\t" +"cantabria"+"\t"+df_gene_c["syn_non_syn_y"][i]+"\n")
        print(df_gene_c["chrom"][i] + "\t" + str(df_gene_c["pos"][i])+ "\t"+ df_gene_c["gene_y"][i] + "\t" +"cantabria"+"\t"+df_gene_c["syn_non_syn_y"][i])
output.close()
