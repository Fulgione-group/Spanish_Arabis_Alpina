import pandas as pd
import statsmodels.formula.api as smf
import numpy as np
import statsmodels.api as sm
df_log=pd.read_csv("output_log_reg_df.tsv",sep="\t")
output=open("statistics_per_gene_polymorphism_g3_g10_filter_low_pvt_var_10_wald0.4184490972257678_gee.tsv","w")
output.write("gene"+"\t"+"OR"+"\t"+"coeff"+"\t"+"pvalues"+"\t"+"bse"+"\t"+"wald_p"+"\n")
genes=df_log["gene"].unique()
R = np.array([[0, 1]])  # select coefficient 1 (second element)
q = np.array([0.4184490972257678])
for gene in genes:
	df_log_gene=df_log[df_log["gene"]==gene]
	if len(df_log_gene[(df_log_gene["set"]=="france") & (df_log_gene["ns"]=="n") ])>=3 and len(df_log_gene[(df_log_gene["set"]=="france") & (df_log_gene["ns"]=="s") ])>=3 and len(df_log_gene[(df_log_gene["set"]=="cantabria") & (df_log_gene["ns"]=="n") ])>=3 and len(df_log_gene[(df_log_gene["set"]=="cantabria") & (df_log_gene["ns"]=="s") ])>=3 and len(df_log_gene[df_log_gene["set"]=="france"])>10 and len(df_log_gene[df_log_gene["set"]=="cantabria"])>10:
		df_log_gene.sort_values(by=['pos'])
		if len(df_log_gene["pos"].value_counts()[df_log_gene["pos"].value_counts()==1])>10:
		        syn=[]
		        for i in df_log_gene["ns"]:
		            if i =="n":
		                syn.append(0)
		            else:
		                syn.append(1)
		        df_log_gene["synonymous"]=syn
		        try:
		            model = sm.GEE.from_formula(
		                "synonymous ~ C(set)",
		                groups="pos",
		                data=df_log_gene,family=sm.families.Binomial()
		            ).fit()
		            wald_result = model.wald_test((R, q), scalar=True)
		            print(wald_result.pvalue)
		            print( gene, np.exp(model.params[1]), model.params[1], model.pvalues[1], model.bse[1])
		            output.write(gene+"\t"+str(np.exp(model.params[1]))+"\t"+str(model.params[1])+"\t"+str(model.pvalues[1])+"\t"+str(model.bse[1])+"\t"+str(wald_result.pvalue)+"\n")
		        except:
		            print("exception")
output.close()
