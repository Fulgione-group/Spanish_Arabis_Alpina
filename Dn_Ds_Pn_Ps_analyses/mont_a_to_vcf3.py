import gzip
import csv

# Load the data file and create a dictionary for quick lookups
with open("gt_data_mont_a_3.txt", "r") as data_file:
    data = data_file.readlines()[1:]
    lookup_dict = {}
    for line in data:
        line = line.rstrip()
        parts = line.split(",")
        key = parts[0] + "," + parts[1]
        lookup_dict[key] = parts[2] + ":"

# Process the VCF file and write the output
with gzip.open("/biodata/dep_coupland/grp_fulgione/common/1000Genomes/GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.b.vcf.gz", "rt") as vcf, \
     open("GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.b.with_mont_a_3.vcf", "w") as new_vcf:

    tsvin = csv.reader(vcf, delimiter="\t")
    c = csv.writer(new_vcf, delimiter="\t", quoting=csv.QUOTE_NONE, escapechar=" ")
    
    for row in tsvin:
        if row[0].startswith("##"):
            c.writerow(row)
        elif row[0] == "#CHROM":
            row.append("mont_a")  # Append the new string before writing
            c.writerow(row)
        else:
            key = row[0] + "," + row[1]
            new_element = lookup_dict[key]+ ":".join(['.'] * (len(row[9].split(":"))-1))
            row.append(new_element)
            c.writerow(row)
