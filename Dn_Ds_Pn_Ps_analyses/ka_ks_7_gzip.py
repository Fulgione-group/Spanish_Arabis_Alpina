import vcf
import csv
import sys
from vcf.parser import Reader, Writer
import gzip
# grep '^#\|missense_variant\|synonymous_variant' annotated_vcf.vcf > mis_syn.txt 

file = sys.argv[1]
outgroup = sys.argv[2]
out = sys.argv[3]
vcf_reader = vcf.Reader(filename=file)
gene_name_index = vcf_reader.infos['ANN'][3].split("|").index(" Gene_Name ")
geneid_index =vcf_reader.infos['ANN'][3].split("|").index(" Gene_ID ")
annotation_index = vcf_reader.infos['ANN'][3].split("|").index(" Annotation ")
feature_id_index = vcf_reader.infos['ANN'][3].split("|").index(" Feature_ID ")
#key = geneid, value = [ka, ks, gene_name]
def add_gene(record):
    output=open(out,"a")
    #output.write("chrom\tpos\tdn\tds\tgene_id\tsm\n")
    #df.columns=["chr","pos","daf_non_syn","daf_syn","gene","ns","ps"]
    curr_gene = str(record.INFO['ANN']).split('|')[geneid_index]
    r=record.INFO['ANN'][0]
    gene = str(r).split('|')[geneid_index]
    if str(r).split('|')[annotation_index] == "synonymous_variant" :
        if record.genotype(outgroup)['GT']=="0/0":
            if record.num_het + 2 * record.num_hom_ref<=2:
                 output.write(record.CHROM+"\t"+str(record.POS)+"\t"+"0"+'\t'+ str(record.num_het + 2 * record.num_hom_alt)+"\t"+str(r).split('|')[gene_name_index]+"\t"+"s"+"\t"+"s\n")
            else:
                output.write(record.CHROM+"\t"+str(record.POS)+"\t"+"0"+'\t'+ str(record.num_het + 2 * record.num_hom_alt)+"\t"+str(r).split('|')[gene_name_index]+"\t"+"s"+"\t"+"p\n")
        if record.genotype(outgroup)['GT']=="1/1":
            if record.num_het + 2 * record.num_hom_alt<=2:
                output.write(record.CHROM+"\t"+str(record.POS)+"\t"+"0"+'\t'+ str(record.num_het + 2 * record.num_hom_ref)+"\t"+str(r).split('|')[gene_name_index]+"\t"+"s"+"\t"+"s\n")
            else:
                output.write(record.CHROM+"\t"+str(record.POS)+"\t"+"0"+'\t'+ str(record.num_het + 2 * record.num_hom_ref)+"\t"+str(r).split('|')[gene_name_index]+"\t"+"s"+"\t"+"p\n")
    if str(r).split('|')[annotation_index] == "missense_variant" :
        if record.genotype(outgroup)['GT']=="0/0":
            if record.num_het + 2 * record.num_hom_ref<=2:
                output.write(record.CHROM+"\t"+str(record.POS)+"\t"+str(record.num_het + 2 * record.num_hom_alt)+"\t"+"0"+"\t"+str(r).split('|')[gene_name_index]+"\t"+"n"+"\t"+"s\n")
            else:
                output.write(record.CHROM+"\t"+str(record.POS)+"\t"+str(record.num_het + 2 * record.num_hom_alt)+"\t"+"0"+"\t"+str(r).split('|')[gene_name_index]+"\t"+"n"+"\tp\n")
        if record.genotype(outgroup)['GT']=="1/1":
            if record.num_het + 2 * record.num_hom_alt<=2:
                output.write(record.CHROM+"\t"+str(record.POS)+"\t"+str(record.num_het + 2 * record.num_hom_ref)+"\t"+"0"+"\t"+str(r).split('|')[gene_name_index]+"\t"+"n"+"\t"+"s\n")
            else:
                output.write(record.CHROM+"\t"+str(record.POS)+"\t"+str(record.num_het + 2 * record.num_hom_ref)+"\t"+"0"+"\t"+str(r).split('|')[gene_name_index]+"\t"+"n"+"\tp\n")
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

##mehak sharma
