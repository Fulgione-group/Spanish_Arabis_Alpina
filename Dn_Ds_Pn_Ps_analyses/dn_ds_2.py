import argparse
import pandas as pd

parser = argparse.ArgumentParser(prog="dn_ds", description="calculates the dnds from the calculated freq")

parser.add_argument('filename')
parser.add_argument("-o", '--output')
args = parser.parse_args()
file=args.filename
out=args.output

df=pd.read_csv(file, sep="\t",header=None)
df.columns=["chr","pos","daf_non_syn","daf_syn","gene","ns","ps"]

genes=[]
output=open(out,"a")
output.write("gene\tDn\tDs\tPn\tPs\n")
for i in df["gene"]:
    if i not in genes:
        dfi=df[df["gene"]==i]
        dfins=dfi[dfi["ns"]=="n"]
        dfins=dfins[dfins["ps"]=="s"]
        dfins=dfins[dfins["daf_non_syn"]!=0] ##this is the missing data
        ns=len(dfins)
        dfinp=dfi[dfi["ns"]=="n"]
        dfinp=dfinp[dfinp["ps"]=="p"]
        dfinp=dfinp[dfinp["daf_non_syn"]!=0]
        np=len(dfinp)
        dfiss=dfi[dfi["ns"]=="s"]
        dfiss=dfiss[dfiss["ps"]=="s"]
        dfiss=dfiss[dfiss["daf_syn"]!=0]
        ss=len(dfiss)
        dfisp=dfi[dfi["ns"]=="s"]
        dfisp=dfisp[dfisp["ps"]=="p"]
        dfisp=dfisp[dfisp["daf_syn"]!=0]
        sp=len(dfisp)
        genes.append(i)
        output=open(out,"a")
        output.write(i+"\t"+ str(ns)+"\t"+str(ss)+"\t"+str(np)+"\t"+str(sp)+"\n")
        print(i+"\t"+ str(ns)+"\t"+str(ss)+"\t"+str(np)+"\t"+str(sp)+"\n")
        output.close()
