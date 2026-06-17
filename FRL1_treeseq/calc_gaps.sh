#!/bin/bash
for chr in {1..8}; do
awk -v chr="chr${chr}" '
NR==1 {next}
NR==2 {prev=$1; next}
{
    realdist=$1-prev
    ratio=$2/realdist

    if (realdist >= 1000 && ratio < 0.01)
        print chr, prev, $1, realdist, $2, ratio

    prev=$1
}
' OFS="\t" relate_popsize_CAN_FR_chr${chr}.dist \
> low_density_gaps_chr${chr}.bed
done
