import vcf
import csv
import sys
from vcf.parser import Reader, Writer

file = sys.argv[1]
outgroup = sys.argv[2]
out=sys.argv[3]
vcf_reader = vcf.Reader(open(file,"r"))
gene_name_index = vcf_reader.infos['ANN'][3].split("|").index(" Gene_Name ")
geneid_index =vcf_reader.infos['ANN'][3].split("|").index(" Gene_ID ")
annotation_index = vcf_reader.infos['ANN'][3].split("|").index(" Annotation ")
feature_id_index = vcf_reader.infos['ANN'][3].split("|").index(" Feature_ID ")

def add_gene(record):
    output=open(out,"a")
    #output.write("chrom\tpos\tpi\tgene_id\n")
    curr_gene = str(record.INFO['ANN']).split('|')[geneid_index]
    r=record.INFO['ANN'][0]
    gene = str(r).split('|')[geneid_index]
    comparisons=0
    differences=0
    if str(r).split('|')[annotation_index] == "synonymous_variant" :
        try:
            output.write(record.CHROM+"\t"+str(record.POS)+"\t"+str(record.nucl_diversity)+"\t"+str(r).split('|')[gene_name_index]+"\t"+"s"+"\n")
        except ZeroDivisionError:
            output.write(record.CHROM+"\t"+str(record.POS)+"\t"+"-"+"\t"+str(r).split('|')[gene_name_index]+"\t"+"s"+"\n")
    if str(r).split('|')[annotation_index] == "missense_variant" :
        try:
            output.write(record.CHROM+"\t"+str(record.POS)+"\t"+str(record.nucl_diversity)+"\t"+str(r).split('|')[gene_name_index]+"\t"+"n"+"\n")
        except ZeroDivisionError:
            output.write(record.CHROM+"\t"+str(record.POS)+"\t"+"-"+"\t"+str(r).split('|')[gene_name_index]+"\t"+"n"+"\n")
    output.close()


heterozygous_in_defined_outgroup=0
multiallelic=0
for record in vcf_reader:
    if record.genotype(outgroup) in record.get_hets():
        heterozygous_in_defined_outgroup=heterozygous_in_defined_outgroup+1
        continue
    elif len(record.alleles[1:])>1:
        multiallelic=multiallelic+1
        continue
    else:
        add_gene(record)
print("heterozygous_in_defined_outgroup",heterozygous_in_defined_outgroup)
print("multiallelic",multiallelic)
