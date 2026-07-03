zcat  GROUP1.vcf.gz | \
awk -v seed=123 '
BEGIN{srand(seed); OFS="\t"}

# meta header
/^##/ {print; next}

# sample header: randomly choose haplotype 1 or 2 per sample
/^#CHROM/ {
    for(i=10; i<=NF; i++) {
        hap[i] = int(rand()*2) + 1
        # print sample choice to stderr
        print $i, "hap" hap[i] > "/dev/stderr"
    }
    print
    next
}

# variant lines
{
    fmt=$9
    for(i=10; i<=NF; i++) {
        split($i, fields, ":")
        gt = fields[1]

        # split phased genotype
        n = split(gt, a, /\|/)

        if(n == 2) {
            fields[1] = a[hap[i]]
        } else {
            fields[1] = "."
        }

        # rebuild sample field
        $i = fields[1]
        for(j=2; j<=length(fields); j++) {
            $i = $i ":" fields[j]
        }
    }
    print
}' | bgzip > GROUP1.haploidized.vcf.gz

tabix -p vcf GROUP1.haploidized.vcf.gz
