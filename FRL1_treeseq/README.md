# Analysis of tree sequences at the FRL1 locus using msprime

This analysis builds on the Relate-derived tree sequences described in the **Demography** directory.

## Identification of low-coverage regions

Regions with poor marker density were identified by flagging intervals where:

- The distance between consecutive SNPs is ≥ 1 kb, and
- Effective coverage is < 1% of the physical distance.

Script:

`calc_gaps.sh`

## Calculation of tree statistics

For all trees, the following statistics were calculated:

- Root age
- Tree span
- Colless index
- Root balance
- SNP density
- Concordance between Cantabrian and French samples at the root split

Script:

`calc_tree_stats_no_lowdensity_with_concordance_notmrca.py`

Additionaly stats
- Root to child ratio

Script:

`calc_tree_stats_no_lowdensity_with_concordance_root_child_ratio.py`

## Identification of long haplotypes separating Cantabrian and French lineages

Consecutive trees were merged into haplotype blocks when they met the following criteria:

- At least 95% agreement between Cantabrian and French samples at the root split
- Minimum root-child size of 152 samples on both sides of the root split
- No overlap with previously identified low-density regions

Script:

`anc_mut_chr_CAN_FR_with_genes.sh`
