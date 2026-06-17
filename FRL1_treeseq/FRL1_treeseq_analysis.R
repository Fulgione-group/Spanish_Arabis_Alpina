########## Compare FRL1 trees to all trees, age, collesIndex, rootBalance, density, 

tree_summaries <- read.table("All_trees_summary_2_no_lowdensity_with_CANFR_concordance.tsv", header=T)
head(tree_summaries)
tree_summaries$treeStart<- floor(tree_summaries$treeStart +1)
tree_summaries$treeEnd <- floor(tree_summaries$treeEnd)
tree_summaries$span <- tree_summaries$treeEnd - tree_summaries$treeStart +1
tree_summaries$timeRoot <-tree_summaries$timeRoot*1.5

tree_summaries$density <- tree_summaries$numMutations/tree_summaries$span
tree_summaries$density


trees_FRL1_haplotype <- tree_summaries[tree_summaries$treeStart>=12608954 & tree_summaries$treeEnd <=12620546 & tree_summaries$chromosome=="chr8",]

minB<-min(trees_FRL1_haplotype$rootBalance)
maxB<-max(trees_FRL1_haplotype$rootBalance)
tree_summaries_balanced <-  tree_summaries[tree_summaries$rootBalance >= minB & tree_summaries$rootBalance <=maxB,]

trees_FRL1_haplotype$collessIndex_quantile <- sapply(
  trees_FRL1_haplotype$collessIndex,
  function(x) sum(tree_summaries$collessIndex <= x) / nrow(tree_summaries)
)
trees_FRL1_haplotype$timeRoot_quantile <- sapply(
  trees_FRL1_haplotype$timeRoot,
  function(x) sum(tree_summaries$timeRoot >= x) / nrow(tree_summaries)
)

trees_FRL1_haplotype$span_quantile <- sapply(
  trees_FRL1_haplotype$span,
  function(x) sum(tree_summaries$span >= x) / nrow(tree_summaries)
)


trees_FRL1_haplotype$rootBalance_quantile <- sapply(
  trees_FRL1_haplotype$rootBalance,
  function(x) sum(tree_summaries$rootBalance >= x) / nrow(tree_summaries)
)

trees_FRL1_haplotype$density_quantile <- sapply(
  trees_FRL1_haplotype$density,
  function(x) sum(tree_summaries$density >= x) / nrow(tree_summaries)
)


### compare among similarly balanced tres
trees_FRL1_haplotype$timeRoot_quantile_balanced <- sapply(
  trees_FRL1_haplotype$timeRoot,
  function(x) sum(tree_summaries_balanced$timeRoot >= x) / nrow(tree_summaries_balanced)
)
hist(tree_summaries_balanced$timeRoot)

trees_FRL1_haplotype$span_quantile_balanced <- sapply(
  trees_FRL1_haplotype$span,
  function(x) sum(tree_summaries_balanced$span >= x) / nrow(tree_summaries_balanced)
)




write.table(trees_FRL1_haplotype, file = "trees_FRL1_haplotypes_compared2all_trees_noLow_density.txt", sep = "\t", col.names=T, row.names = F, quote = F)


############### Analyze the consecutive trees that seperate haplotypes 

balanced_trees <- read.table("balanced_rootsplit_runs_95agreement_no_lowdensity_with_root_and_branch_ages.tsv", sep="\t", header = T)

median(balanced_trees$mean_root_child_branch, n=100)
hist(balanced_trees$mean_root_child_branch, n=100)

FRL1_trees<-balanced_trees[balanced_trees$start==12608954,]

sum(balanced_trees$span_bp > FRL1_trees$span_bp)/nrow(balanced_trees)

