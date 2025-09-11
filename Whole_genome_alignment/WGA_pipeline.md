# Whole Genome Alignment Pipeline  
*A. alpina* vs. *A. montbretiana*

This document describes the workflow used to produce syntenic positions of *A. montbretiana* projected onto *A. alpina* coordinates.  
It contains both the **pre-RepeatFiller alignment** and the **post-RepeatFiller filtering and conversion** steps.

---

## Pre-RepeatFiller Steps

1. **Chunk alignment with lastz**  
Split the query genome into chunks (e.g. 10 Mb). Align each chunk:

```
lastz --gap=400,30 --gappedthresh=3000 --seed=12of19 --inner=2200 --format=axt+   Arabis_alpina.MPIPZ.version_5.1.chr.all.2bit[multiple]   Amontbretiana_genomic_chr1_0${i}.2bit > Amontbretiana_genomic_chr1_0${i}.axt
```

2. **Convert AXT to chains (per chunk)**  
```
axtChain -linearGap=medium   Amontbretiana_genomic_chr1_0${i}.axt   Arabis_alpina.MPIPZ.version_5.1.chr.all.2bit   Amontbretiana_genomic_chr1_0${i}.2bit   Amontbretiana_genomic_chr1_0${i}.chain
```

3. **Merge and sort all chains**  
```
chainMergeSort -inputList=all_chains.txt -tempDir=/scratch > all_merged_sorted.chain
```

4. **Patch chains (sensitive realignments)**  
```
patchChain.perl all_merged_sorted.chain   Arabis_alpina.MPIPZ.version_5.1.chr.all.2bit   GCA_001484125.2_MPIPZ_Amontbretiana_3.1_genomic.2bit   Arabis_alpina.MPIPZ.version_5.1.chr.all.genome.sizes   GCA_001484125.2_MPIPZ_Amontbretiana_3.1_genomic.genome.sizes   -lastzParameters "--gap=400,30 --gappedthresh=3000 --seed=12of19 --inner=2200 --format=axt+"   -outputDir ./patchChain -numJobs 100
```

5. **RepeatFiller**  
```
RepeatFiller.py -c all_merged_sorted.chain   -T2 Arabis_alpina.MPIPZ.version_5.1.chr.all.2bit   -Q2 GCA_001484125.2_MPIPZ_Amontbretiana_3.1_genomic.2bit   --lastzParameters 'gap=400,30 gappedthresh=3000 seed=12of19 inner=2200 format=axt+'   --workdir /scratch   --output all_withPatch_repeatFilled_merged_sorted.chain
```

---

## Post-RepeatFiller Steps

1. **ChainCleaner (optional)**  
```
chainCleaner all_withPatch_repeatFilled_merged_sorted.chain   Arabis_alpina.MPIPZ.version_5.1.chr.all.2bit   GCA_001484125.2_MPIPZ_Amontbretiana_3.1_genomic.2bit   all_merged_sorted_patched_repfilled_cleaned.chain   all_merged_sorted_patched_repfilled_cleaned.removed   -tSizes=Arabis_alpina.MPIPZ.version_5.1.chr.all.genome.sizes   -qSizes=GCA_001484125.2_MPIPZ_Amontbretiana_3.1_genomic.genome.sizes   -linearGap=medium
```

2. **Merge and sort patched chains (no IDs)**  
```
chainMergeSort -inputList=all_chains_with_newAlignment.txt -tempDir=/scratch > all_withPatch_merged_sorted_noID.chain
```

3. **Prepare prenet**  
```
chainPreNet all_withPatch_merged_sorted_noID.chain   Arabis_alpina.MPIPZ.version_5.1.chr.all.genome.sizes   GCA_001484125.2_MPIPZ_Amontbretiana_3.1_genomic.genome.sizes   all_withPatch_merged_sorted_noID.chain.prenet
```

4. **Netting**  
```
chainNet -rescore -linearGap=medium   -tNibDir=Arabis_alpina.MPIPZ.version_5.1.chr.all.2bit   -qNibDir=GCA_001484125.2_MPIPZ_Amontbretiana_3.1_genomic.2bit   all_withPatch_merged_sorted_noID.chain.prenet   Arabis_alpina.MPIPZ.version_5.1.chr.all.genome.sizes   GCA_001484125.2_MPIPZ_Amontbretiana_3.1_genomic.genome.sizes   all_withPatch_merged_sorted_alpina_noID.net   all_withPatch_merged_sorted_mont_noID.net
```

5. **Keep only syntenic nets**  
```
netSyntenic all_withPatch_merged_sorted_alpina_noID.net all_withPatch_merged_sorted_alpina_noID_syntenic.net
netSyntenic all_withPatch_merged_sorted_mont_noID.net   all_withPatch_merged_sorted_mont_noID_syntenic.net
```

6. **Convert nets to AXT**  
```
netToAxt all_withPatch_merged_sorted_alpina_noID_syntenic.net   all_withPatch_merged_sorted_noID.chain.prenet   Arabis_alpina.MPIPZ.version_5.1.chr.all.2bit   GCA_001484125.2_MPIPZ_Amontbretiana_3.1_genomic.2bit   all_withPatch_merged_noID_sorted_alpina_syntenic.axt
```

7. **Convert AXT to MAF**  
```
axtToMaf all_withPatch_merged_noID_sorted_alpina_syntenic.axt   Arabis_alpina.MPIPZ.version_5.1.chr.all.genome.sizes   GCA_001484125.2_MPIPZ_Amontbretiana_3.1_genomic.genome.sizes   all_withPatch_merged_noID_sorted_alpina_syntenic.maf
```

8. **Convert MAF to BED (example: chr1)**  
```
maf_to_bed.py   -i all_withPatch_merged_noID_sorted_alpina_syntenic_montID_alpinaID_chr1.maf.gz   -r alpina -c chr1 | sort -k1,1 -k2,2n   > all_withPatch_merged_noID_sorted_alpina_syntenic_montID_alpinaID_chr1.bed
```

9. **Optional: build FASTA from BED**  
```
awk ... # your make_fasta_from_bed.sh script logic
```
