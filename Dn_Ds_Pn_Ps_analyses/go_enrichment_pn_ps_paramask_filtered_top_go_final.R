library(topGO)
BPterms <- ls(GOBPTerm) # loading the biological processes terms from topgo
library(tidyr)
find("annFUN.db")
annFUN.db
data(geneList)
gene2go_df <- read.table("GO_alpina_extended_homology_based_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gene2go_df)=c("GeneID","GOterms")
geneID2GO <- setNames(strsplit(gene2go_df$GOterms, " "), gene2go_df$GeneID)
head(geneID2GO)

geneList <- read.table("significant_genes_pn_ps_paramask.txt", header = FALSE, stringsAsFactors = FALSE)[,1]
allGenes <- unique(gene2go_df$GeneID)
find("geneList")

geneSelection <- factor(as.integer(allGenes %in% geneList))
names(geneSelection) <- allGenes

geneSelection
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneSelection,
              geneSel = function(p) p == 1,  # Function to select significant genes
              description = "GO analysis",
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)
find("annFUN.gene2GO")
GOdata


resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")


resultFisher.elim
resultClassic

allRes <- GenTable(GOdata, classicFisher = resultClassic,
                   elimFisher = resultFisher.elim,
                   orderBy = "elimFisher", ranksOf = "classicFisher", topNodes = 10)
allRes

go_genes <- genesInTerm(GOdata, "GO:0006468")
go_genes <- genesInTerm(GOdata, "GO:0010583")
geneSelection[names(geneSelection)%in%go_genes$`GO:0006468` & geneSelection=="1"]
geneSelection[names(geneSelection)%in%go_genes$`GO:0010583` & geneSelection=="1"]
go_genes <- genesInTerm(GOdata, "GO:0120010")
geneSelection[names(geneSelection)%in%go_genes$`GO:0120010` & geneSelection=="1"]
go_genes <- genesInTerm(GOdata, "GO:0007178")
geneSelection[names(geneSelection)%in%go_genes$`GO:0007178` & geneSelection=="1"]
go_genes <- genesInTerm(GOdata, "GO:0045806")
geneSelection[names(geneSelection)%in%go_genes$`GO:0045806` & geneSelection=="1"]
geneSelection[geneSelection=="1"]
