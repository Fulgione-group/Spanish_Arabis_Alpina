import vcf
import csv
import sys
from vcf.parser import Reader, Writer 
import pandas as pd

sample_to_add=open("all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_chrall_sorted_mlines_noGaps_uppercases.fasta","r")
sample=sample_to_add.readlines()
sample.remove("\n")
ref=open("../annotation/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta","r")
ref=ref.readlines()
data=open("gt_data_mont_a_3.txt","a")
data.write("chr,pos,gt,ref,sample\n")
data.close()
cordinates=[]
for i in ref:
    if i.startswith(">"):
        cordinates.append(ref.index(i))
print(cordinates[0:8])
cordinates_sample=[]
for i in sample:
    if i.startswith(">"):
        cordinates_sample.append(sample.index(i))
print(cordinates_sample[0:8])

for c in [0,1,2,3,4,5,6,7]:
    info=ref[cordinates[c]+1:cordinates[c+1]]
    for i in range(len(info)):
        info[i]=info[i].rstrip()
    if c!=7:
        sample_i=sample[cordinates_sample[c]+1:cordinates_sample[c+1]]
    if c==7:
        sample_i=sample[cordinates_sample[c]+1:]
    for i in range(len(sample_i)):
        sample_i[i]=sample_i[i].rstrip()
    sample_i="".join(sample_i)
    info="".join(info)
    print(len(sample_i),len(info))
    for i in range(min(len(info),len(sample_i))):
        if info[i]==sample_i[i]:
            data=open("gt_data_mont_a_3.txt","a")
            data.write("chr"+str(c+1)+","+str(i+1)+","+"0/0"+","+info[i]+","+sample_i[i]+"\n")
            data.close()
        else:
            if sample_i[i]=="N":
                data=open("gt_data_mont_a_3.txt","a")
                data.write("chr"+str(c+1)+","+str(i+1)+","+"./."+","+info[i]+","+sample_i[i]+"\n")
                data.close()
            else:
                data=open("gt_data_mont_a_3.txt","a")
                data.write("chr"+str(c+1)+","+str(i+1)+","+"1/1"+","+info[i]+","+sample_i[i]+"\n")
                data.close()

