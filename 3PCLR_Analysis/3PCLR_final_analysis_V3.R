library(qqman)
library(ggplot2)
library(cowplot)
library(grid)
library(Rgraphviz)
library(topGO)
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

###loadMask
mask<- read.table(file = "ParaMask_Cantabria_and_French_30PerMiss_v3_EM2.7.1_EMresults.chrall.finalClass.multicopy.ZeroBasedStart.SNPable.3col.sorted.merged.bed", sep = "\t")
head(mask)
colnames(mask)<- c("chr", "start", "end")
mask$start<- mask$start+1


mask_FR<- read.table(file = "ParaMask_French_clean_EM_2.7.1_EMresults.chrall.finalClass.bed", sep = "\t", header = T)
head(mask_FR)
mask_CN<- read.table(file = "Cantabria_30PerMiss_v3_EM2.7.1_EMresults.chrall.finalClass.bed", sep = "\t", header = T)



#3.1
Set3.1_xpclr_chrall<- read.table(file = "3PCLR_v3/3PCLR_G1_G2_G3/Set3.1_xpclr_chrall_new.txt", sep = "\t", header=T)

c1<- sort(Set3.1_xpclr_chrall$X3PCLR.Anc.std, decreasing = T)[nrow(Set3.1_xpclr_chrall)*0.001]

#3.2

Set3.2_xpclr_chrall<- read.table(file = "3PCLR_v3/3PCLR_G1_ES17_G4//Set3.2_xpclr_chrall_new.txt", sep = "\t", header=T)

c1<- sort(Set3.2_xpclr_chrall$X3PCLR.Anc.std, decreasing = T)[nrow(Set3.2_xpclr_chrall)*0.001]

#3.3

Set3.3_xpclr_chrall<- read.table(file = "3PCLR_v3/3PCLR_G2_ES17_G4/Set3.3_xpclr_chrall_new.txt", sep = "\t", header=T)

c1<- sort(Set3.3_xpclr_chrall$X3PCLR.Anc.std, decreasing = T)[nrow(Set3.3_xpclr_chrall)*0.001]


# #2.1
#
Set2.1_xpclr_chrall<- read.table(file = "3PCLR_v3/3PCLR_G1_G2_FR/Set2.1_xpclr_chrall_new.txt", sep = "\t", header=T)
#
# c1<-sort(Set2.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.001*length(Set2.1_xpclr_chrall$X3PCLR.Anc.std)]
#
#2.2

Set2.2_xpclr_chrall<- read.table(file = "3PCLR_v3/3PCLR_G1_ES17_FR/Set2.2_xpclr_chrall_new.txt", sep = "\t", header=T)

#
# c1<-sort(Set2.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.001*length(Set2.2_xpclr_chrall$X3PCLR.Anc.std)]
#
# #2.3
Set2.3_xpclr_chrall<- read.table(file = "3PCLR_v3/3PCLR_G2_ES17_FR/Set2.3_xpclr_chrall_new.txt", sep = "\t", header=T)
#
# c1<-sort(Set2.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.001*length(Set2.3_xpclr_chrall$X3PCLR.Anc.std)]

#merge sets for plotting
Set2.1_xpclr_chrall$set<-2.1
Set2.2_xpclr_chrall$set<-2.2
Set2.3_xpclr_chrall$set<-2.3
Set3.1_xpclr_chrall$set<-3.1
Set3.2_xpclr_chrall$set<-3.2
Set3.3_xpclr_chrall$set<-3.3
Set3.3_xpclr_chrall
SetAll_xpclr_chrall <- rbind(Set2.1_xpclr_chrall, Set2.2_xpclr_chrall, Set2.3_xpclr_chrall,Set3.1_xpclr_chrall ,Set3.2_xpclr_chrall ,Set3.3_xpclr_chrall)
SetALL_xpclr_chrall <- rbind(Set2.1_xpclr_chrall, Set2.2_xpclr_chrall, Set2.3_xpclr_chrall,Set3.1_xpclr_chrall ,Set3.2_xpclr_chrall ,Set3.3_xpclr_chrall)

###find overlapping genes


Set2.1_xpclr_chrall$start<- Set2.1_xpclr_chrall$PhysStart
Set2.1_xpclr_chrall$end<- Set2.1_xpclr_chrall$PhysEnd
Set2.1_xpclr_chrall$end[1] <-floor(Set2.1_xpclr_chrall$PhysPos[1]+ (Set2.1_xpclr_chrall$PhysPos[2]-Set2.1_xpclr_chrall$PhysPos[1])/2)
for(i in 2:nrow(Set2.1_xpclr_chrall)){
  Set2.1_xpclr_chrall$start[i] <- ceiling(0.5+Set2.1_xpclr_chrall$PhysPos[i]+ (Set2.1_xpclr_chrall$PhysPos[i-1]-Set2.1_xpclr_chrall$PhysPos[i])/2)
  Set2.1_xpclr_chrall$end[i] <- floor(Set2.1_xpclr_chrall$PhysPos[i]+ (Set2.1_xpclr_chrall$PhysPos[i+1]-Set2.1_xpclr_chrall$PhysPos[i])/2)
}
Set2.1_xpclr_chrall$end[nrow(Set2.1_xpclr_chrall)] <- Set2.1_xpclr_chrall$PhysEnd[nrow(Set2.1_xpclr_chrall)]

Set2.1_xpclr_chrall$start[Set2.1_xpclr_chrall$Chr=="chr2"][1] <-Set2.1_xpclr_chrall$PhysStart[Set2.1_xpclr_chrall$Chr=="chr2"][1]
Set2.1_xpclr_chrall$start[Set2.1_xpclr_chrall$Chr=="chr3"][1] <-Set2.1_xpclr_chrall$PhysStart[Set2.1_xpclr_chrall$Chr=="chr3"][1]
Set2.1_xpclr_chrall$start[Set2.1_xpclr_chrall$Chr=="chr4"][1] <-Set2.1_xpclr_chrall$PhysStart[Set2.1_xpclr_chrall$Chr=="chr4"][1]
Set2.1_xpclr_chrall$start[Set2.1_xpclr_chrall$Chr=="chr5"][1] <-Set2.1_xpclr_chrall$PhysStart[Set2.1_xpclr_chrall$Chr=="chr5"][1]
Set2.1_xpclr_chrall$start[Set2.1_xpclr_chrall$Chr=="chr6"][1] <-Set2.1_xpclr_chrall$PhysStart[Set2.1_xpclr_chrall$Chr=="chr6"][1]
Set2.1_xpclr_chrall$start[Set2.1_xpclr_chrall$Chr=="chr7"][1] <-Set2.1_xpclr_chrall$PhysStart[Set2.1_xpclr_chrall$Chr=="chr7"][1]
Set2.1_xpclr_chrall$start[Set2.1_xpclr_chrall$Chr=="chr8"][1] <-Set2.1_xpclr_chrall$PhysStart[Set2.1_xpclr_chrall$Chr=="chr8"][1]
Set2.1_xpclr_chrall$end[Set2.1_xpclr_chrall$Chr=="chr1"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr1",])] <- Set2.1_xpclr_chrall$PhysEnd[Set2.1_xpclr_chrall$Chr=="chr1"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr1",])]
Set2.1_xpclr_chrall$end[Set2.1_xpclr_chrall$Chr=="chr2"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr2",])] <- Set2.1_xpclr_chrall$PhysEnd[Set2.1_xpclr_chrall$Chr=="chr2"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr2",])]
Set2.1_xpclr_chrall$end[Set2.1_xpclr_chrall$Chr=="chr3"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr3",])] <- Set2.1_xpclr_chrall$PhysEnd[Set2.1_xpclr_chrall$Chr=="chr3"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr3",])]
Set2.1_xpclr_chrall$end[Set2.1_xpclr_chrall$Chr=="chr4"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr4",])] <- Set2.1_xpclr_chrall$PhysEnd[Set2.1_xpclr_chrall$Chr=="chr4"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr4",])]
Set2.1_xpclr_chrall$end[Set2.1_xpclr_chrall$Chr=="chr5"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr5",])] <- Set2.1_xpclr_chrall$PhysEnd[Set2.1_xpclr_chrall$Chr=="chr5"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr5",])]
Set2.1_xpclr_chrall$end[Set2.1_xpclr_chrall$Chr=="chr6"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr6",])] <- Set2.1_xpclr_chrall$PhysEnd[Set2.1_xpclr_chrall$Chr=="chr6"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr6",])]
Set2.1_xpclr_chrall$end[Set2.1_xpclr_chrall$Chr=="chr7"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr7",])] <- Set2.1_xpclr_chrall$PhysEnd[Set2.1_xpclr_chrall$Chr=="chr7"][nrow(Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr7",])]




Set2.2_xpclr_chrall$start<- Set2.2_xpclr_chrall$PhysStart
Set2.2_xpclr_chrall$end<- Set2.2_xpclr_chrall$PhysEnd
Set2.2_xpclr_chrall$end[1] <-floor(Set2.2_xpclr_chrall$PhysPos[1]+ (Set2.2_xpclr_chrall$PhysPos[2]-Set2.2_xpclr_chrall$PhysPos[1])/2)
for(i in 2:nrow(Set2.2_xpclr_chrall)){
  Set2.2_xpclr_chrall$start[i] <- ceiling(0.5+Set2.2_xpclr_chrall$PhysPos[i]+ (Set2.2_xpclr_chrall$PhysPos[i-1]-Set2.2_xpclr_chrall$PhysPos[i])/2)
  Set2.2_xpclr_chrall$end[i] <- floor(Set2.2_xpclr_chrall$PhysPos[i]+ (Set2.2_xpclr_chrall$PhysPos[i+1]-Set2.2_xpclr_chrall$PhysPos[i])/2)
}
Set2.2_xpclr_chrall$start[Set2.2_xpclr_chrall$Chr=="chr2"][1] <-Set2.2_xpclr_chrall$PhysStart[Set2.2_xpclr_chrall$Chr=="chr2"][1]
Set2.2_xpclr_chrall$start[Set2.2_xpclr_chrall$Chr=="chr3"][1] <-Set2.2_xpclr_chrall$PhysStart[Set2.2_xpclr_chrall$Chr=="chr3"][1]
Set2.2_xpclr_chrall$start[Set2.2_xpclr_chrall$Chr=="chr4"][1] <-Set2.2_xpclr_chrall$PhysStart[Set2.2_xpclr_chrall$Chr=="chr4"][1]
Set2.2_xpclr_chrall$start[Set2.2_xpclr_chrall$Chr=="chr5"][1] <-Set2.2_xpclr_chrall$PhysStart[Set2.2_xpclr_chrall$Chr=="chr5"][1]
Set2.2_xpclr_chrall$start[Set2.2_xpclr_chrall$Chr=="chr6"][1] <-Set2.2_xpclr_chrall$PhysStart[Set2.2_xpclr_chrall$Chr=="chr6"][1]
Set2.2_xpclr_chrall$start[Set2.2_xpclr_chrall$Chr=="chr7"][1] <-Set2.2_xpclr_chrall$PhysStart[Set2.2_xpclr_chrall$Chr=="chr7"][1]
Set2.2_xpclr_chrall$start[Set2.2_xpclr_chrall$Chr=="chr8"][1] <-Set2.2_xpclr_chrall$PhysStart[Set2.2_xpclr_chrall$Chr=="chr8"][1]
Set2.2_xpclr_chrall$end[Set2.2_xpclr_chrall$Chr=="chr1"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr1",])] <- Set2.2_xpclr_chrall$PhysEnd[Set2.2_xpclr_chrall$Chr=="chr1"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr1",])]
Set2.2_xpclr_chrall$end[Set2.2_xpclr_chrall$Chr=="chr2"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr2",])] <- Set2.2_xpclr_chrall$PhysEnd[Set2.2_xpclr_chrall$Chr=="chr2"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr2",])]
Set2.2_xpclr_chrall$end[Set2.2_xpclr_chrall$Chr=="chr3"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr3",])] <- Set2.2_xpclr_chrall$PhysEnd[Set2.2_xpclr_chrall$Chr=="chr3"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr3",])]
Set2.2_xpclr_chrall$end[Set2.2_xpclr_chrall$Chr=="chr4"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr4",])] <- Set2.2_xpclr_chrall$PhysEnd[Set2.2_xpclr_chrall$Chr=="chr4"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr4",])]
Set2.2_xpclr_chrall$end[Set2.2_xpclr_chrall$Chr=="chr5"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr5",])] <- Set2.2_xpclr_chrall$PhysEnd[Set2.2_xpclr_chrall$Chr=="chr5"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr5",])]
Set2.2_xpclr_chrall$end[Set2.2_xpclr_chrall$Chr=="chr6"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr6",])] <- Set2.2_xpclr_chrall$PhysEnd[Set2.2_xpclr_chrall$Chr=="chr6"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr6",])]
Set2.2_xpclr_chrall$end[Set2.2_xpclr_chrall$Chr=="chr7"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr7",])] <- Set2.2_xpclr_chrall$PhysEnd[Set2.2_xpclr_chrall$Chr=="chr7"][nrow(Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr7",])]





Set2.3_xpclr_chrall$start<- Set2.3_xpclr_chrall$PhysStart
Set2.3_xpclr_chrall$end<- Set2.3_xpclr_chrall$PhysEnd
Set2.3_xpclr_chrall$end[1] <-floor(Set2.3_xpclr_chrall$PhysPos[1]+ (Set2.3_xpclr_chrall$PhysPos[2]-Set2.3_xpclr_chrall$PhysPos[1])/2)
for(i in 2:nrow(Set2.3_xpclr_chrall)){
  Set2.3_xpclr_chrall$start[i] <- ceiling(0.5+Set2.3_xpclr_chrall$PhysPos[i]+ (Set2.3_xpclr_chrall$PhysPos[i-1]-Set2.3_xpclr_chrall$PhysPos[i])/2)
  Set2.3_xpclr_chrall$end[i] <- floor(Set2.3_xpclr_chrall$PhysPos[i]+ (Set2.3_xpclr_chrall$PhysPos[i+1]-Set2.3_xpclr_chrall$PhysPos[i])/2)
}
Set2.3_xpclr_chrall$end[nrow(Set2.3_xpclr_chrall)] <- Set2.3_xpclr_chrall$PhysEnd[nrow(Set2.3_xpclr_chrall)]
Set2.3_xpclr_chrall$start[Set2.3_xpclr_chrall$Chr=="chr2"][1] <-Set2.3_xpclr_chrall$PhysStart[Set2.3_xpclr_chrall$Chr=="chr2"][1]
Set2.3_xpclr_chrall$start[Set2.3_xpclr_chrall$Chr=="chr3"][1] <-Set2.3_xpclr_chrall$PhysStart[Set2.3_xpclr_chrall$Chr=="chr3"][1]
Set2.3_xpclr_chrall$start[Set2.3_xpclr_chrall$Chr=="chr4"][1] <-Set2.3_xpclr_chrall$PhysStart[Set2.3_xpclr_chrall$Chr=="chr4"][1]
Set2.3_xpclr_chrall$start[Set2.3_xpclr_chrall$Chr=="chr5"][1] <-Set2.3_xpclr_chrall$PhysStart[Set2.3_xpclr_chrall$Chr=="chr5"][1]
Set2.3_xpclr_chrall$start[Set2.3_xpclr_chrall$Chr=="chr6"][1] <-Set2.3_xpclr_chrall$PhysStart[Set2.3_xpclr_chrall$Chr=="chr6"][1]
Set2.3_xpclr_chrall$start[Set2.3_xpclr_chrall$Chr=="chr7"][1] <-Set2.3_xpclr_chrall$PhysStart[Set2.3_xpclr_chrall$Chr=="chr7"][1]
Set2.3_xpclr_chrall$start[Set2.3_xpclr_chrall$Chr=="chr8"][1] <-Set2.3_xpclr_chrall$PhysStart[Set2.3_xpclr_chrall$Chr=="chr8"][1]
Set2.3_xpclr_chrall$end[Set2.3_xpclr_chrall$Chr=="chr1"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr1",])] <- Set2.3_xpclr_chrall$PhysEnd[Set2.3_xpclr_chrall$Chr=="chr1"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr1",])]
Set2.3_xpclr_chrall$end[Set2.3_xpclr_chrall$Chr=="chr2"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr2",])] <- Set2.3_xpclr_chrall$PhysEnd[Set2.3_xpclr_chrall$Chr=="chr2"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr2",])]
Set2.3_xpclr_chrall$end[Set2.3_xpclr_chrall$Chr=="chr3"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr3",])] <- Set2.3_xpclr_chrall$PhysEnd[Set2.3_xpclr_chrall$Chr=="chr3"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr3",])]
Set2.3_xpclr_chrall$end[Set2.3_xpclr_chrall$Chr=="chr4"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr4",])] <- Set2.3_xpclr_chrall$PhysEnd[Set2.3_xpclr_chrall$Chr=="chr4"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr4",])]
Set2.3_xpclr_chrall$end[Set2.3_xpclr_chrall$Chr=="chr5"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr5",])] <- Set2.3_xpclr_chrall$PhysEnd[Set2.3_xpclr_chrall$Chr=="chr5"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr5",])]
Set2.3_xpclr_chrall$end[Set2.3_xpclr_chrall$Chr=="chr6"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr6",])] <- Set2.3_xpclr_chrall$PhysEnd[Set2.3_xpclr_chrall$Chr=="chr6"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr6",])]
Set2.3_xpclr_chrall$end[Set2.3_xpclr_chrall$Chr=="chr7"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr7",])] <- Set2.3_xpclr_chrall$PhysEnd[Set2.3_xpclr_chrall$Chr=="chr7"][nrow(Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr7",])]




Set3.1_xpclr_chrall$start<- Set3.1_xpclr_chrall$PhysStart
Set3.1_xpclr_chrall$end<- Set3.1_xpclr_chrall$PhysEnd
Set3.1_xpclr_chrall$end[1] <-floor(Set3.1_xpclr_chrall$PhysPos[1]+ (Set3.1_xpclr_chrall$PhysPos[2]-Set3.1_xpclr_chrall$PhysPos[1])/2)
for(i in 2:nrow(Set3.1_xpclr_chrall)){
  Set3.1_xpclr_chrall$start[i] <- ceiling(0.5+Set3.1_xpclr_chrall$PhysPos[i]+ (Set3.1_xpclr_chrall$PhysPos[i-1]-Set3.1_xpclr_chrall$PhysPos[i])/2)
  Set3.1_xpclr_chrall$end[i] <- floor(Set3.1_xpclr_chrall$PhysPos[i]+ (Set3.1_xpclr_chrall$PhysPos[i+1]-Set3.1_xpclr_chrall$PhysPos[i])/2)
}
Set3.1_xpclr_chrall$end[nrow(Set3.1_xpclr_chrall)] <- Set3.1_xpclr_chrall$PhysEnd[nrow(Set3.1_xpclr_chrall)]
Set3.1_xpclr_chrall$start[Set3.1_xpclr_chrall$Chr=="chr2"][1] <-Set3.1_xpclr_chrall$PhysStart[Set3.1_xpclr_chrall$Chr=="chr2"][1]
Set3.1_xpclr_chrall$start[Set3.1_xpclr_chrall$Chr=="chr3"][1] <-Set3.1_xpclr_chrall$PhysStart[Set3.1_xpclr_chrall$Chr=="chr3"][1]
Set3.1_xpclr_chrall$start[Set3.1_xpclr_chrall$Chr=="chr4"][1] <-Set3.1_xpclr_chrall$PhysStart[Set3.1_xpclr_chrall$Chr=="chr4"][1]
Set3.1_xpclr_chrall$start[Set3.1_xpclr_chrall$Chr=="chr5"][1] <-Set3.1_xpclr_chrall$PhysStart[Set3.1_xpclr_chrall$Chr=="chr5"][1]
Set3.1_xpclr_chrall$start[Set3.1_xpclr_chrall$Chr=="chr6"][1] <-Set3.1_xpclr_chrall$PhysStart[Set3.1_xpclr_chrall$Chr=="chr6"][1]
Set3.1_xpclr_chrall$start[Set3.1_xpclr_chrall$Chr=="chr7"][1] <-Set3.1_xpclr_chrall$PhysStart[Set3.1_xpclr_chrall$Chr=="chr7"][1]
Set3.1_xpclr_chrall$start[Set3.1_xpclr_chrall$Chr=="chr8"][1] <-Set3.1_xpclr_chrall$PhysStart[Set3.1_xpclr_chrall$Chr=="chr8"][1]
Set3.1_xpclr_chrall$end[Set3.1_xpclr_chrall$Chr=="chr1"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr1",])] <- Set3.1_xpclr_chrall$PhysEnd[Set3.1_xpclr_chrall$Chr=="chr1"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr1",])]
Set3.1_xpclr_chrall$end[Set3.1_xpclr_chrall$Chr=="chr2"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr2",])] <- Set3.1_xpclr_chrall$PhysEnd[Set3.1_xpclr_chrall$Chr=="chr2"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr2",])]
Set3.1_xpclr_chrall$end[Set3.1_xpclr_chrall$Chr=="chr3"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3",])] <- Set3.1_xpclr_chrall$PhysEnd[Set3.1_xpclr_chrall$Chr=="chr3"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3",])]
Set3.1_xpclr_chrall$end[Set3.1_xpclr_chrall$Chr=="chr4"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr4",])] <- Set3.1_xpclr_chrall$PhysEnd[Set3.1_xpclr_chrall$Chr=="chr4"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr4",])]
Set3.1_xpclr_chrall$end[Set3.1_xpclr_chrall$Chr=="chr5"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr5",])] <- Set3.1_xpclr_chrall$PhysEnd[Set3.1_xpclr_chrall$Chr=="chr5"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr5",])]
Set3.1_xpclr_chrall$end[Set3.1_xpclr_chrall$Chr=="chr6"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr6",])] <- Set3.1_xpclr_chrall$PhysEnd[Set3.1_xpclr_chrall$Chr=="chr6"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr6",])]
Set3.1_xpclr_chrall$end[Set3.1_xpclr_chrall$Chr=="chr7"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr7",])] <- Set3.1_xpclr_chrall$PhysEnd[Set3.1_xpclr_chrall$Chr=="chr7"][nrow(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr7",])]


Set3.2_xpclr_chrall$start<- Set3.2_xpclr_chrall$PhysStart
Set3.2_xpclr_chrall$end<- Set3.2_xpclr_chrall$PhysEnd
Set3.2_xpclr_chrall$end[1] <-floor(Set3.2_xpclr_chrall$PhysPos[1]+ (Set3.2_xpclr_chrall$PhysPos[2]-Set3.2_xpclr_chrall$PhysPos[1])/2)
for(i in 2:nrow(Set3.2_xpclr_chrall)){
  Set3.2_xpclr_chrall$start[i] <- ceiling(0.5+Set3.2_xpclr_chrall$PhysPos[i]+ (Set3.2_xpclr_chrall$PhysPos[i-1]-Set3.2_xpclr_chrall$PhysPos[i])/2)
  Set3.2_xpclr_chrall$end[i] <- floor(Set3.2_xpclr_chrall$PhysPos[i]+ (Set3.2_xpclr_chrall$PhysPos[i+1]-Set3.2_xpclr_chrall$PhysPos[i])/2)
}
Set3.2_xpclr_chrall$end[nrow(Set3.2_xpclr_chrall)] <- Set3.2_xpclr_chrall$PhysEnd[nrow(Set3.2_xpclr_chrall)]
Set3.2_xpclr_chrall$start[Set3.2_xpclr_chrall$Chr=="chr2"][1] <-Set3.2_xpclr_chrall$PhysStart[Set3.2_xpclr_chrall$Chr=="chr2"][1]
Set3.2_xpclr_chrall$start[Set3.2_xpclr_chrall$Chr=="chr3"][1] <-Set3.2_xpclr_chrall$PhysStart[Set3.2_xpclr_chrall$Chr=="chr3"][1]
Set3.2_xpclr_chrall$start[Set3.2_xpclr_chrall$Chr=="chr4"][1] <-Set3.2_xpclr_chrall$PhysStart[Set3.2_xpclr_chrall$Chr=="chr4"][1]
Set3.2_xpclr_chrall$start[Set3.2_xpclr_chrall$Chr=="chr5"][1] <-Set3.2_xpclr_chrall$PhysStart[Set3.2_xpclr_chrall$Chr=="chr5"][1]
Set3.2_xpclr_chrall$start[Set3.2_xpclr_chrall$Chr=="chr6"][1] <-Set3.2_xpclr_chrall$PhysStart[Set3.2_xpclr_chrall$Chr=="chr6"][1]
Set3.2_xpclr_chrall$start[Set3.2_xpclr_chrall$Chr=="chr7"][1] <-Set3.2_xpclr_chrall$PhysStart[Set3.2_xpclr_chrall$Chr=="chr7"][1]
Set3.2_xpclr_chrall$start[Set3.2_xpclr_chrall$Chr=="chr8"][1] <-Set3.2_xpclr_chrall$PhysStart[Set3.2_xpclr_chrall$Chr=="chr8"][1]
Set3.2_xpclr_chrall$end[Set3.2_xpclr_chrall$Chr=="chr1"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr1",])] <- Set3.2_xpclr_chrall$PhysEnd[Set3.2_xpclr_chrall$Chr=="chr1"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr1",])]
Set3.2_xpclr_chrall$end[Set3.2_xpclr_chrall$Chr=="chr2"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr2",])] <- Set3.2_xpclr_chrall$PhysEnd[Set3.2_xpclr_chrall$Chr=="chr2"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr2",])]
Set3.2_xpclr_chrall$end[Set3.2_xpclr_chrall$Chr=="chr3"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr3",])] <- Set3.2_xpclr_chrall$PhysEnd[Set3.2_xpclr_chrall$Chr=="chr3"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr3",])]
Set3.2_xpclr_chrall$end[Set3.2_xpclr_chrall$Chr=="chr4"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr4",])] <- Set3.2_xpclr_chrall$PhysEnd[Set3.2_xpclr_chrall$Chr=="chr4"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr4",])]
Set3.2_xpclr_chrall$end[Set3.2_xpclr_chrall$Chr=="chr5"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr5",])] <- Set3.2_xpclr_chrall$PhysEnd[Set3.2_xpclr_chrall$Chr=="chr5"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr5",])]
Set3.2_xpclr_chrall$end[Set3.2_xpclr_chrall$Chr=="chr6"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr6",])] <- Set3.2_xpclr_chrall$PhysEnd[Set3.2_xpclr_chrall$Chr=="chr6"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr6",])]
Set3.2_xpclr_chrall$end[Set3.2_xpclr_chrall$Chr=="chr7"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr7",])] <- Set3.2_xpclr_chrall$PhysEnd[Set3.2_xpclr_chrall$Chr=="chr7"][nrow(Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr7",])]



Set3.3_xpclr_chrall$start<- Set3.3_xpclr_chrall$PhysStart
Set3.3_xpclr_chrall$end<- Set3.3_xpclr_chrall$PhysEnd
Set3.3_xpclr_chrall$end[1] <-floor(Set3.3_xpclr_chrall$PhysPos[1]+ (Set3.3_xpclr_chrall$PhysPos[2]-Set3.3_xpclr_chrall$PhysPos[1])/2)
for(i in 2:nrow(Set3.3_xpclr_chrall)){
  Set3.3_xpclr_chrall$start[i] <- ceiling(0.5+Set3.3_xpclr_chrall$PhysPos[i]+ (Set3.3_xpclr_chrall$PhysPos[i-1]-Set3.3_xpclr_chrall$PhysPos[i])/2)
  Set3.3_xpclr_chrall$end[i] <- floor(Set3.3_xpclr_chrall$PhysPos[i]+ (Set3.3_xpclr_chrall$PhysPos[i+1]-Set3.3_xpclr_chrall$PhysPos[i])/2)
}
Set3.3_xpclr_chrall$end[nrow(Set3.3_xpclr_chrall)] <- Set3.3_xpclr_chrall$PhysEnd[nrow(Set3.3_xpclr_chrall)]
Set3.3_xpclr_chrall$start[Set3.3_xpclr_chrall$Chr=="chr2"][1] <-Set3.3_xpclr_chrall$PhysStart[Set3.3_xpclr_chrall$Chr=="chr2"][1]
Set3.3_xpclr_chrall$start[Set3.3_xpclr_chrall$Chr=="chr3"][1] <-Set3.3_xpclr_chrall$PhysStart[Set3.3_xpclr_chrall$Chr=="chr3"][1]
Set3.3_xpclr_chrall$start[Set3.3_xpclr_chrall$Chr=="chr4"][1] <-Set3.3_xpclr_chrall$PhysStart[Set3.3_xpclr_chrall$Chr=="chr4"][1]
Set3.3_xpclr_chrall$start[Set3.3_xpclr_chrall$Chr=="chr5"][1] <-Set3.3_xpclr_chrall$PhysStart[Set3.3_xpclr_chrall$Chr=="chr5"][1]
Set3.3_xpclr_chrall$start[Set3.3_xpclr_chrall$Chr=="chr6"][1] <-Set3.3_xpclr_chrall$PhysStart[Set3.3_xpclr_chrall$Chr=="chr6"][1]
Set3.3_xpclr_chrall$start[Set3.3_xpclr_chrall$Chr=="chr7"][1] <-Set3.3_xpclr_chrall$PhysStart[Set3.3_xpclr_chrall$Chr=="chr7"][1]
Set3.3_xpclr_chrall$start[Set3.3_xpclr_chrall$Chr=="chr8"][1] <-Set3.3_xpclr_chrall$PhysStart[Set3.3_xpclr_chrall$Chr=="chr8"][1]
Set3.3_xpclr_chrall$end[Set3.3_xpclr_chrall$Chr=="chr1"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr1",])] <- Set3.3_xpclr_chrall$PhysEnd[Set3.3_xpclr_chrall$Chr=="chr1"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr1",])]
Set3.3_xpclr_chrall$end[Set3.3_xpclr_chrall$Chr=="chr2"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr2",])] <- Set3.3_xpclr_chrall$PhysEnd[Set3.3_xpclr_chrall$Chr=="chr2"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr2",])]
Set3.3_xpclr_chrall$end[Set3.3_xpclr_chrall$Chr=="chr3"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr3",])] <- Set3.3_xpclr_chrall$PhysEnd[Set3.3_xpclr_chrall$Chr=="chr3"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr3",])]
Set3.3_xpclr_chrall$end[Set3.3_xpclr_chrall$Chr=="chr4"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr4",])] <- Set3.3_xpclr_chrall$PhysEnd[Set3.3_xpclr_chrall$Chr=="chr4"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr4",])]
Set3.3_xpclr_chrall$end[Set3.3_xpclr_chrall$Chr=="chr5"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr5",])] <- Set3.3_xpclr_chrall$PhysEnd[Set3.3_xpclr_chrall$Chr=="chr5"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr5",])]
Set3.3_xpclr_chrall$end[Set3.3_xpclr_chrall$Chr=="chr6"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr6",])] <- Set3.3_xpclr_chrall$PhysEnd[Set3.3_xpclr_chrall$Chr=="chr6"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr6",])]
Set3.3_xpclr_chrall$end[Set3.3_xpclr_chrall$Chr=="chr7"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr7",])] <- Set3.3_xpclr_chrall$PhysEnd[Set3.3_xpclr_chrall$Chr=="chr7"][nrow(Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr7",])]







### select genes for different branches for GO


c<-sort(Set3.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.002*length(Set3.1_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set3.1_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.0005*length(Set3.1_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.00001*length(Set3.1_xpclr_chrall$X3PCLR.Anc.std)]
# c <- 14.5
gff<- read.table("Arabis_alpina_mpipz_v5.1_annotation.genes.bed", header = F, sep="\t")
head(gff)
Set3.1_xpclr_chrall_Anc_select<- Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$X3PCLR.Anc.std>=c,]
nrow(Set3.1_xpclr_chrall_Anc_select)
head(gff)
head(Set3.1_xpclr_chrall_Anc_select)
gff$score<- 0
i <-2
gff_select <- c()

for(i in 1:nrow(Set3.1_xpclr_chrall_Anc_select)){
  start <- Set3.1_xpclr_chrall_Anc_select$start[i]
  end <- Set3.1_xpclr_chrall_Anc_select$end[i]
  chr <- Set3.1_xpclr_chrall_Anc_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set3.1_xpclr_chrall_Anc_select$X3PCLR.Anc.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}

gff_select<-gff_select[order(gff_select$score, decreasing = T),]
unique(gff_select$V9)
gff_select
AA3G31360
gff_select[gff_select$V9=="AA3G31360",]
gff_select[gff_select$V9=="AA1G29480",]
gff_select[gff_select$V9=="AA2G06360",]
#MUCI70 paralog
gff_select[gff_select$V9=="AA5G02960",]


gff_select[order(gff_select$score, decreasing = T),]



###selec
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.001*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.002*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.0005*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.00001*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]
c<-sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]

# c <- 14.5
gff<- read.table("Arabis_alpina_mpipz_v5.1_annotation.genes.bed", header = F, sep="\t")
head(gff)
Set3.1_xpclr_chrall_A_select<- Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$X3PCLR.A.std>=c,]
nrow(Set3.1_xpclr_chrall_A_select)
head(gff)
head(Set3.1_xpclr_chrall_A_select)
gff$score<- 0
i <-2
gff_select <- c()

for(i in 1:nrow(Set3.1_xpclr_chrall_A_select)){
  start <- Set3.1_xpclr_chrall_A_select$start[i]
  end <- Set3.1_xpclr_chrall_A_select$end[i]
  chr <- Set3.1_xpclr_chrall_A_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set3.1_xpclr_chrall_A_select$X3PCLR.A.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
gff_select

gff_select[gff_select$V9=="AA3G31360",]

#VIP4
gff_select[gff_select$V9=="AA8G47450",]
gff_select$score[gff_select$V9=="AA8G47450"]

sum(sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA8G47450"])/length(Set3.1_xpclr_chrall$X3PCLR.A.std)

##

# c<-sort(Set3.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.001*length(Set3.1_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.002*length(Set3.1_xpclr_chrall$X3PCLR.B.std)]
c<-sort(Set3.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set3.1_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.0005*length(Set3.1_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.00001*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]
# c <- 14.5
gff<- read.table("Arabis_alpina_mpipz_v5.1_annotation.genes.bed", header = F, sep="\t")
head(gff)
Set3.1_xpclr_chrall_B_select<- Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$X3PCLR.B.std>=c,]
nrow(Set3.1_xpclr_chrall_B_select)
head(gff)
head(Set3.1_xpclr_chrall_B_select)
gff$score<- 0
i <-2
gff_select<-gff_select[order(gff_select$score, decreasing = T),]
gff_select
# # gff_select <- c()

for(i in 1:nrow(Set3.1_xpclr_chrall_B_select)){
    start <- Set3.1_xpclr_chrall_B_select$start[i]
  end <- Set3.1_xpclr_chrall_B_select$end[i]
  chr <- Set3.1_xpclr_chrall_B_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set3.1_xpclr_chrall_B_select$X3PCLR.B.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
gff_select
gff_select[order(gff_select$score, decreasing = T),]
gff_select[gff_select$V9=="AA6G27810",]
gff_select[gff_select$V9=="AA3G31360",]


sum(sort(Set3.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA1G35160"])/length(Set3.1_xpclr_chrall$X3PCLR.B.std)



##3.2
###

c<-sort(Set3.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set3.2_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.0005*length(Set3.2_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.00001*length(Set3.2_xpclr_chrall$X3PCLR.Anc.std)]
# c <- 14.5
gff<- read.table("Arabis_alpina_mpipz_v5.1_annotation.genes.bed", header = F, sep="\t")
head(gff)
Set3.2_xpclr_chrall_Anc_select<- Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$X3PCLR.Anc.std>=c,]
nrow(Set3.2_xpclr_chrall_Anc_select)
head(gff)
head(Set3.2_xpclr_chrall_Anc_select)
gff$score<- 0
i <-2
# gff_select <- c()

for(i in 1:nrow(Set3.2_xpclr_chrall_Anc_select)){
  start <- Set3.2_xpclr_chrall_Anc_select$start[i]
  end <- Set3.2_xpclr_chrall_Anc_select$end[i]
  chr <- Set3.2_xpclr_chrall_Anc_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set3.2_xpclr_chrall_Anc_select$X3PCLR.Anc.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}

gff_select<-gff_select[order(gff_select$score, decreasing = T),]
unique(gff_select$V9)

gff_select[gff_select$V9=="AA1G29480",]
gff_select[gff_select$V9=="AA3G31360",]

#EZA 1
sum(sort(Set3.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA6G21780"])/length(Set3.2_xpclr_chrall$X3PCLR.Anc.std)
#VIL1
sum(sort(Set3.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA3G25750"])/length(Set3.2_xpclr_chrall$X3PCLR.Anc.std)



###selec
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.001*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.002*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]
c<-sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.0005*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.00001*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]
# c <- 14.5
gff<- read.table("Arabis_alpina_mpipz_v5.1_annotation.genes.bed", header = F, sep="\t")
head(gff)
Set3.2_xpclr_chrall_A_select<- Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$X3PCLR.A.std>=c,]
nrow(Set3.2_xpclr_chrall_A_select)
head(gff)
head(Set3.2_xpclr_chrall_A_select)
gff$score<- 0
i <-2
# # gff_select <- c()

for(i in 1:nrow(Set3.2_xpclr_chrall_A_select)){
  start <- Set3.2_xpclr_chrall_A_select$start[i]
  end <- Set3.2_xpclr_chrall_A_select$end[i]
  chr <- Set3.2_xpclr_chrall_A_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set3.2_xpclr_chrall_A_select$X3PCLR.A.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
gff_select
# gff_select[gff_select$V9=="AA6G27810",]
gff_select[gff_select$V9=="AA6G27810",]
gff_select[gff_select$V9=="AA6G27830",]
gff_select[gff_select$V9=="AA3G31360",]

#EZA 1
sum(sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA6G21780"])/length(Set3.2_xpclr_chrall$X3PCLR.A.std)
#VIL1
sum(sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA3G25750"])/length(Set3.2_xpclr_chrall$X3PCLR.A.std)

# c<-sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.002*length(Set3.2_xpclr_chrall$X3PCLR.B.std)]
c<-sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set3.2_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.0005*length(Set3.2_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.00001*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]
# c <- 14.5
gff<- read.table("Arabis_alpina_mpipz_v5.1_annotation.genes.bed", header = F, sep="\t")
head(gff)
Set3.2_xpclr_chrall_B_select<- Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$X3PCLR.B.std>=c,]
nrow(Set3.2_xpclr_chrall_B_select)
head(gff)
head(Set3.2_xpclr_chrall_B_select)
gff$score<- 0
i <-2
gff_select<-gff_select[order(gff_select$score, decreasing = T),]
gff_select
# # gff_select <- c()

for(i in 1:nrow(Set3.2_xpclr_chrall_B_select)){
  start <- Set3.2_xpclr_chrall_B_select$start[i]
  end <- Set3.2_xpclr_chrall_B_select$end[i]
  chr <- Set3.2_xpclr_chrall_B_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set3.2_xpclr_chrall_B_select$X3PCLR.B.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
gff_select
gff_select[order(gff_select$score, decreasing = T),]
gff_select[gff_select$V9=="AA6G27810",]
gff_select[gff_select$V9=="AA3G31360",]

  #EZA 1
sum(sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA6G21780"])/length(Set3.2_xpclr_chrall$X3PCLR.B.std)
#VIL1
sum(sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA3G25750"])/length(Set3.2_xpclr_chrall$X3PCLR.B.std)
sum(sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA1G35160"])/length(Set3.2_xpclr_chrall$X3PCLR.B.std)

##3.3
###

c<-sort(Set3.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set3.3_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set3.3_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.0005*length(Set3.3_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.00001*length(Set3.3_xpclr_chrall$X3PCLR.Anc.std)]
# c <- 14.5
gff<- read.table("Arabis_alpina_mpipz_v5.1_annotation.genes.bed", header = F, sep="\t")
head(gff)
Set3.3_xpclr_chrall_Anc_select<- Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$X3PCLR.Anc.std>=c,]
nrow(Set3.3_xpclr_chrall_Anc_select)
head(gff)
head(Set3.3_xpclr_chrall_Anc_select)
gff$score<- 0
i <-2
# # gff_select <- c()

for(i in 1:nrow(Set3.3_xpclr_chrall_Anc_select)){
  start <- Set3.3_xpclr_chrall_Anc_select$start[i]
  end <- Set3.3_xpclr_chrall_Anc_select$end[i]
  chr <- Set3.3_xpclr_chrall_Anc_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set3.3_xpclr_chrall_Anc_select$X3PCLR.Anc.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}

gff_select<-gff_select[order(gff_select$score, decreasing = T),]
unique(gff_select$V9)
gff_select[gff_select$V9=="AA1G29480",]
gff_select[gff_select$V9=="AA3G31360",]



#EZA 1
sum(sort(Set3.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA6G21780"])/length(Set3.3_xpclr_chrall$X3PCLR.Anc.std)
#VIL1
sum(sort(Set3.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA3G25750"])/length(Set3.3_xpclr_chrall$X3PCLR.Anc.std)

###selec
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.002*length(Set3.3_xpclr_chrall$X3PCLR.A.std)]
c<-sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set3.3_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.002*length(Set3.3_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set3.3_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.0005*length(Set3.3_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.00001*length(Set3.3_xpclr_chrall$X3PCLR.A.std)]
# c <- 14.5
gff <- read.table("Arabis_alpina_mpipz_v5.1_annotation.genes.bed", header = F, sep="\t")
head(gff)
Set3.3_xpclr_chrall_A_select<- Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$X3PCLR.A.std>=c,]
nrow(Set3.3_xpclr_chrall_A_select)
head(gff)
head(Set3.3_xpclr_chrall_A_select)
gff$score<- 0
i <-2
# # gff_select <- c()

for(i in 1:nrow(Set3.3_xpclr_chrall_A_select)){
  start <- Set3.3_xpclr_chrall_A_select$start[i]
  end <- Set3.3_xpclr_chrall_A_select$end[i]
  chr <- Set3.3_xpclr_chrall_A_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set3.3_xpclr_chrall_A_select$X3PCLR.A.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
gff_select
gff_select[gff_select$V9=="AA6G27810",]
gff_select[gff_select$V9=="AA6G27830",]
gff_select[gff_select$V9=="AA3G31360",]

#EZA 1
sum(sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA6G21780"])/length(Set3.3_xpclr_chrall$X3PCLR.A.std)
#VIL1
sum(sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA3G25750"])/length(Set3.3_xpclr_chrall$X3PCLR.A.std)

c<-sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[(0.005*length(Set3.3_xpclr_chrall$X3PCLR.B.std))]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[(0.002*length(Set3.3_xpclr_chrall$X3PCLR.B.std))]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set3.3_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.0005*length(Set3.3_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.00001*length(Set3.3_xpclr_chrall$X3PCLR.A.std)]
# c <- 14.5
gff<- read.table("Arabis_alpina_mpipz_v5.1_annotation.genes.bed", header = F, sep="\t")
head(gff)
Set3.3_xpclr_chrall_B_select<- Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$X3PCLR.B.std>=c,]
nrow(Set3.3_xpclr_chrall_B_select)
head(gff)
head(Set3.3_xpclr_chrall_B_select)
gff$score<- 0
i <-2
gff_select<-gff_select[order(gff_select$score, decreasing = T),]
gff_select
# # gff_select <- c()

for(i in 1:nrow(Set3.3_xpclr_chrall_B_select)){
  start <- Set3.3_xpclr_chrall_B_select$start[i]
  end <- Set3.3_xpclr_chrall_B_select$end[i]
  chr <- Set3.3_xpclr_chrall_B_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set3.3_xpclr_chrall_B_select$X3PCLR.B.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
gff_select
gff_select[order(gff_select$score, decreasing = T),]
gff_select[gff_select$V9=="AA6G27810",]
gff_select[gff_select$V9=="AA3G31360",]


sum(sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA6G21780"])/length(Set3.3_xpclr_chrall$X3PCLR.B.std)
#VIL1
sum(sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)>=gff_select$score[gff_select$V9=="AA3G25750"])/length(Set3.3_xpclr_chrall$X3PCLR.B.std)




##
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.001*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.002*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]
c<-sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.0005*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.00001*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]
# c <- 14.5
gff<- read.table("Arabis_alpina_mpipz_v5.1_annotation.genes.bed", header = F, sep="\t")
head(gff)
Set2.1_xpclr_chrall_A_select<- Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$X3PCLR.A.std>=c,]
nrow(Set2.1_xpclr_chrall_A_select)
head(gff)
head(Set2.1_xpclr_chrall_A_select)
gff$score<- 0
i <-2
gff_select <- c()

for(i in 1:nrow(Set2.1_xpclr_chrall_A_select)){
  start <- Set2.1_xpclr_chrall_A_select$start[i]
  end <- Set2.1_xpclr_chrall_A_select$end[i]
  chr <- Set2.1_xpclr_chrall_A_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set2.1_xpclr_chrall_A_select$X3PCLR.A.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}

#top hit GUF1
nrow(gff_select)
gff_select[gff_select$V9=="AA3G31360",]

gff_select[order(gff_select$score, decreasing = T),]
#


# c<-sort(Set2.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.002*length(Set2.1_xpclr_chrall$X3PCLR.B.std)]
c<-sort(Set2.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.001*length(Set2.1_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set2.1_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.0005*length(Set2.1_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.00001*length(Set2.1_xpclr_chrall$X3PCLR.B.std)]

Set2.1_xpclr_chrall_B_select<- Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$X3PCLR.B.std>=c,]
nrow(Set2.1_xpclr_chrall_A_select)
head(gff)
head(Set2.1_xpclr_chrall_A_select)
gff$score<- 0
i <-2
# # gff_select <- c()

for(i in 1:nrow(Set2.1_xpclr_chrall_B_select)){
  start <- Set2.1_xpclr_chrall_B_select$start[i]
  end <- Set2.1_xpclr_chrall_B_select$end[i]
  chr <- Set2.1_xpclr_chrall_B_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set2.1_xpclr_chrall_B_select$X3PCLR.B.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
#top hits IFRD1, ACA1
nrow(gff_select)
gff_select[order(gff_select$score, decreasing = T),]







## private branches set2.2 and 2.3

# c<-sort(Set2.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.002*length(Set2.2_xpclr_chrall$X3PCLR.A.std)]
c<-sort(Set2.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set2.2_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.001*length(Set2.2_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set2.2_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.0005*length(Set2.2_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.0001*length(Set2.2_xpclr_chrall$X3PCLR.A.std)]

Set2.2_xpclr_chrall_A_select<- Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$X3PCLR.A.std>=c,]
gff$score<- 0
i <-2
# # gff_select <- c()

for(i in 1:nrow(Set2.2_xpclr_chrall_A_select)){
  start <- Set2.2_xpclr_chrall_A_select$start[i]
  end <- Set2.2_xpclr_chrall_A_select$end[i]
  chr <- Set2.2_xpclr_chrall_A_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set2.2_xpclr_chrall_A_select$X3PCLR.A.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}

nrow(gff_select)
unique(gff_select$V9)

# c<-sort(Set2.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.002*length(Set2.2_xpclr_chrall$X3PCLR.B.std)]
c<-sort(Set2.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set2.2_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.001*length(Set2.2_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set2.2_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.0005*length(Set2.2_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.00001*length(Set2.2_xpclr_chrall$X3PCLR.B.std)]
# c <-14.5
Set2.2_xpclr_chrall_B_select<- Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$X3PCLR.B.std>=c,]
gff$score<- 0
i <-2
# # gff_select <- c()


for(i in 1:nrow(Set2.2_xpclr_chrall_B_select)){
  start <- Set2.2_xpclr_chrall_B_select$start[i]
  end <- Set2.2_xpclr_chrall_B_select$end[i]
  chr <- Set2.2_xpclr_chrall_B_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set2.2_xpclr_chrall_B_select$X3PCLR.B.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
nrow(gff_select)

#top hits KAS2, TAN, PCMP-E83, APP2, DMXL1, TSET, DExH-box ATP-dependent RNA helicase
length(unique(gff_select$V9))

# c<-sort(Set2.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.002*length(Set2.3_xpclr_chrall$X3PCLR.A.std)]
c<-sort(Set2.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set2.3_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.001*length(Set2.3_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set2.3_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.0005*length(Set2.3_xpclr_chrall$X3PCLR.A.std)]
# c<-sort(Set2.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.0001*length(Set2.3_xpclr_chrall$X3PCLR.A.std)]

Set2.3_xpclr_chrall_A_select<- Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$X3PCLR.A.std>=c,]
gff$score<- 0
i <-2
# # gff_select <- c()


for(i in 1:nrow(Set2.3_xpclr_chrall_A_select)){
  start <- Set2.3_xpclr_chrall_A_select$start[i]
  end <- Set2.3_xpclr_chrall_A_select$end[i]
  chr <- Set2.3_xpclr_chrall_A_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set2.3_xpclr_chrall_A_select$X3PCLR.A.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
##SDH7A,
nrow(gff_select)
unique(gff_select$V9)
# c<-sort(Set2.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.002*length(Set2.3_xpclr_chrall$X3PCLR.B.std)]
c<-sort(Set2.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set2.3_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.001*length(Set2.3_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set2.3_xpclr_chrall$X3PCLR.B.std)]
# c<-sort(Set2.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.0005*length(Set2.3_xpclr_chrall$X3PCLR.B.std)]

Set2.3_xpclr_chrall_B_select<- Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$X3PCLR.B.std>=c,]
gff$score<- 0
i <-2
# # gff_select <- c()

for(i in 1:nrow(Set2.3_xpclr_chrall_B_select)){
  start <- Set2.3_xpclr_chrall_B_select$start[i]
  end <- Set2.3_xpclr_chrall_B_select$end[i]
  chr <- Set2.3_xpclr_chrall_B_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set2.3_xpclr_chrall_B_select$X3PCLR.B.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
nrow(gff_select)
unique(gff_select$V9)
sort(gff_select$score)
gff_select[order(gff_select$score, decreasing = T),]

##anc



##anc set 2

# c<-sort(Set2.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.002*length(Set2.1_xpclr_chrall$X3PCLR.Anc.std)]
c<-sort(Set2.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set2.1_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set2.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.0005*length(Set2.1_xpclr_chrall$X3PCLR.Anc.std)]

Set2.1_xpclr_chrall_Anc_select<- Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$X3PCLR.Anc.std>=c,]
gff_select <- c()


for(i in 1:nrow(Set2.1_xpclr_chrall_Anc_select)){
  start <- Set2.1_xpclr_chrall_Anc_select$start[i]
  end <- Set2.1_xpclr_chrall_Anc_select$end[i]
  chr <- Set2.1_xpclr_chrall_Anc_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set2.1_xpclr_chrall_Anc_select$X3PCLR.Anc.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
nrow(gff_select)

# c<-sort(Set2.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.002*length(Set2.2_xpclr_chrall$X3PCLR.Anc.std)]
c<-sort(Set2.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set2.2_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set2.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set2.2_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set2.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.0005*length(Set2.2_xpclr_chrall$X3PCLR.Anc.std)]

Set2.2_xpclr_chrall_Anc_select<- Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$X3PCLR.Anc.std>=c,]
gff$score<- 0
i <-2
# # gff_select <- c()

for(i in 1:nrow(Set2.2_xpclr_chrall_Anc_select)){
  start <- Set2.2_xpclr_chrall_Anc_select$start[i]
  end <- Set2.2_xpclr_chrall_Anc_select$end[i]
  chr <- Set2.2_xpclr_chrall_Anc_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set2.2_xpclr_chrall_Anc_select$X3PCLR.Anc.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}

nrow(gff_select)


# c<-sort(Set2.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.002*length(Set2.3_xpclr_chrall$X3PCLR.Anc.std)]
c<-sort(Set2.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set2.3_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set2.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set2.3_xpclr_chrall$X3PCLR.Anc.std)]
# c<-sort(Set2.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.0005*length(Set2.3_xpclr_chrall$X3PCLR.Anc.std)]

# c<-sort(Set2.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.001*length(Set2.1_xpclr_chrall$X3PCLR.Anc.std)]

Set2.3_xpclr_chrall_Anc_select<- Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$X3PCLR.Anc.std>=c,]


for(i in 1:nrow(Set2.3_xpclr_chrall_Anc_select)){
  start <- Set2.3_xpclr_chrall_Anc_select$start[i]
  end <- Set2.3_xpclr_chrall_Anc_select$end[i]
  chr <- Set2.3_xpclr_chrall_Anc_select$Chr[i]
  genes<- gff[gff$V1== chr & gff$V4 <= end & gff$V5 >=start,]
  if(nrow(genes)>=1){
    genes$score <- Set2.3_xpclr_chrall_Anc_select$X3PCLR.Anc.std[i]
    gff_select <- rbind(gff_select, genes)
  }
}
nrow(gff_select)
gff_select[order(gff_select$score, decreasing = T),]


### GO enrichment with TopGO
gene2GO_new<- readMappings(file = "GO_alpina_extended_homology_based_genes.txt", sep = "\t",IDsep = " ")
head(gene2GO_new)

genes<- read.table("GO_alpina_extended_homology_based_genes.txt", sep="\t",header=F)

all_genes<-rep(0, length(genes[,1]))
names(all_genes)<- genes[,1]
head(all_genes)


select_genes<- gff_select[, c(9,10)]
# select_genes<- Stat_tab_CAN[Stat_tab_CAN$Pn_Ps_FR_fisher_pval<=0.05,][,c(1,11)]
# select_genes <- as.data.frame(cbind(sig_P$gene,sig_P$Pn_Ps_FR_fisher_pval))
# select_genes$V2 <- as.numeric(select_genes$V2)
# nrow(Stat_tab_CAN[Stat_tab_CAN$Pn_Ps_FR_fisher_pval<=0.05,][,c(1,11)])


head(select_genes)
colnames(select_genes) <- c("V1", "V2")
# select_genes[, 1] <- gsub("T", "G", select_genes[, 1])

length(all_genes[names(all_genes)%in%select_genes[,1]])

select_genes_uniq<- c()
for(i in unique(select_genes[,1])){
  select_genes_uniq<- rbind(select_genes_uniq, c(i,mean(select_genes[select_genes$V1==i,2])))
}
select_genes_uniq<- as.data.frame(select_genes_uniq)
select_genes_uniq$V2<- as.numeric(select_genes_uniq$V2)

nrow(select_genes_uniq)
head(select_genes_uniq)
select_genes_uniq[select_genes_uniq$V1=="AA3G32640",]
sum(select_genes_uniq[,1]%in%names(all_genes))

length(all_genes[names(all_genes) %in% select_genes_uniq[,1]])
# all_genes[names(all_genes) %in% select_genes_uniq[,1]] <- select_genes_uniq[,2][select_genes_uniq[,1]%in%names(all_genes)]
tmp<-select_genes_uniq[select_genes_uniq[,1]%in%names(all_genes),]
tmp2 <- tmp$V2
names(tmp2) <- tmp$V1
tmp2 <- tmp2[names(all_genes[names(all_genes) %in% select_genes_uniq[,1]])]

all_genes[names(all_genes) %in% select_genes_uniq[,1]] <- tmp2


# #
# gene2GO_new<- gene2GO_new[-grep(x=names(gene2GO_new),pattern ="AAs")]

# Your list of GO terms per gene

# GO term to check
# go_term_to_check <- "GO:0007268"
# go_term_to_check <- "GO:0046928"
# Check which genes contain the GO term
# contains_go <- sapply(gene2GO_new, function(x) go_term_to_check %in% x)

# Show the result
# sum(contains_go)


tailXPCLR <- function(allScore) {
  return(allScore > 0)
}

#159 for 2.23, 109 for 2.1,
length(names(all_genes[all_genes>0]))
# names(all_genes)<- genes[,1]
GOdata_xpclr <- new("topGOdata", ontology = "BP", allGenes = all_genes,annot = annFUN.gene2GO, gene2GO =gene2GO_new ,geneSel=tailXPCLR)


resultFisher.elim <- runTest(GOdata_xpclr, algorithm = "elim", statistic = "fisher")
resultClassic <- runTest(GOdata_xpclr, algorithm = "classic", statistic = "fisher")




par(cex=0.15)
showSigOfNodes(GOdata_xpclr, score(resultFisher.elim), firstSigNodes = 5, useInfo ='def',.NO.CHAR = 200, plotFunction = GOplot)

gene_to_GO <- genesInTerm(GOdata_xpclr)


allRes_elim <- GenTable(GOdata_xpclr,   elimFisher = resultFisher.elim, orderBy = "elimFisher", ranksOf = "elimFisher", topNodes = 1000)
head(allRes_elim)


# filter for either minimum of 2 or 3 genes for each category
allRes_select <- allRes_elim[allRes_elim$Significant> 1,]
allRes_elim[allRes_elim$Significant> 2,]
# allRes_select_classic <- allRes_elim[allRes_classic$Significant> 1,]
allRes_select
allRes_select$ratio <- allRes_select$Significant/allRes_select$Annotated
allRes_select[ 1:40,]
allRes_select$elimFisher <- as.numeric(allRes_select$elimFisher)
allRes_select$Significant <- as.numeric(allRes_select$Significant)
allRes_select$ratio <- as.numeric(allRes_select$ratio)


gff_select[order(gff_select$score, decreasing = T),]

sum(gff_select$V9=="AA6G26170")
  # write.table(allRes_select, file = "3P_CLR_v2/All_private_branches_3_GO.txt", quote = F, sep = "\t")
go_genes <- genesInTerm(GOdata_xpclr, "GO:0009414")
select_genes[select_genes$V1%in%go_genes$`GO:0009414`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0009408")
select_genes[select_genes$V1%in%go_genes$`GO:0009408`,]



go_genes <- genesInTerm(GOdata_xpclr, "GO:0071281")
select_genes[select_genes$V1%in%go_genes$`GO:0007291`,]


go_genes <- genesInTerm(GOdata_xpclr, "GO:0070072")
select_genes[select_genes$V1%in%go_genes$`GO:0070072`,]


go_genes <- genesInTerm(GOdata_xpclr, "GO:0009737")
select_genes[select_genes$V1%in%go_genes$`GO:0009737`,]

  go_genes <- genesInTerm(GOdata_xpclr, "GO:0006952")
select_genes[select_genes$V1%in%go_genes$`GO:0006952`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0009409")
select_genes[select_genes$V1%in%go_genes$`GO:0009409`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0055114")
select_genes[select_genes$V1%in%go_genes$`GO:0055114`,]


go_genes <- genesInTerm(GOdata_xpclr, "GO:0006979")
select_genes[select_genes$V1%in%go_genes$`GO:0006979`,]

###check flowering/vern genes

gff_select[order(gff_select$score, decreasing = T),]


allRes_elim[grep(x = allRes_elim$Term,pattern = "flower"),]
allRes_elim[grep(x = allRes_elim$Term,pattern = "vern"),]
allRes_elim[grep(x = allRes_elim$Term,pattern = "flor"),]
go_genes <- genesInTerm(GOdata_xpclr, "GO:0010229")
select_genes[select_genes$V1%in%go_genes$`GO:0010229`,]
allRes_elim[grep(x = allRes_elim$Term,pattern = "vern"),]
allRes_elim[grep(x = allRes_elim$Term,pattern = "salinity"),]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0017144")
select_genes[select_genes$V1%in%go_genes$`GO:0017144`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0048573")
select_genes[select_genes$V1%in%go_genes$`GO:0048573`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0048575")
select_genes[select_genes$V1%in%go_genes$`GO:0048575`,]


go_genes <- genesInTerm(GOdata_xpclr, "GO:2000028")
select_genes[select_genes$V1%in%go_genes$`GO:2000028`,]
go_genes <- genesInTerm(GOdata_xpclr, "GO:0048573")
select_genes[select_genes$V1%in%go_genes$`GO:0048573`,]


go_genes <- genesInTerm(GOdata_xpclr, "GO:0009908")
select_genes[select_genes$V1%in%go_genes$`GO:0009908`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0042538")
select_genes[select_genes$V1%in%go_genes$`GO:0042538`,]
#vern
go_genes <- genesInTerm(GOdata_xpclr, "GO:0010048")
select_genes[select_genes$V1%in%go_genes$`GO:0010048`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0010219")
select_genes[select_genes$V1%in%go_genes$`GO:0010219`,]


# Outgroup FR
go_genes <- genesInTerm(GOdata_xpclr, "GO:0046246")
select_genes[select_genes$V1%in%go_genes$`GO:0046246`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0007291")
select_genes[select_genes$V1%in%go_genes$`GO:0007291`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:1900034")
select_genes[select_genes$V1%in%go_genes$`GO:1900034`,]


allRes_select[ 1:40,]
go_genes <- genesInTerm(GOdata_xpclr, "GO:0045944")
s<-select_genes[select_genes$V1%in%go_genes$`GO:0045944`,]
unique(s$V1)

#anc of W-CAN CE-CAN outgroup
GO.ID                                        Term Annotated Significant Expected elimFisher       ratio
1   GO:0071281               cellular response to iron ion        64           5     0.31   0.000015 0.078125000
2   GO:0071452         cellular response to singlet oxygen         5           2     0.02   0.000230 0.400000000
4   GO:0010358                                leaf shaping        12           2     0.06   0.001500 0.166666667
5   GO:0031670               cellular response to nutrient        12           2     0.06   0.001500 0.166666667
6   GO:0010075               regulation of meristem growth        59           3     0.29   0.002980 0.050847458
7   GO:0046283 anthocyanin-containing compound metaboli...        62           3     0.30   0.003430 0.048387097
#PB CE-CAN outgroup
GO.ID                                        Term Annotated Significant Expected elimFisher      ratio
1   GO:1902476            chloride transmembrane transport         8           3     0.10    0.00012 0.37500000
2   GO:0035652       clathrin-coated vesicle cargo loading         2           2     0.03    0.00017 1.00000000
4   GO:0006582                   melanin metabolic process         3           2     0.04    0.00050 0.66666667
5   GO:0045921           positive regulation of exocytosis         3           2     0.04    0.00050 0.66666667
6   GO:0001676     long-chain fatty acid metabolic process        32           4     0.42    0.00075 0.12500000
7   GO:0043162 ubiquitin-dependent protein catabolic pr...        15           3     0.19    0.00088 0.20000000
8   GO:0050848    regulation of calcium-mediated signaling         4           2     0.05    0.00099 0.50000000
10  GO:0000381 regulation of alternative mRNA splicing,...        19           3     0.25    0.00179 0.15789474
12  GO:0051209 release of sequestered calcium ion into ...         7           2     0.09    0.00337 0.28571429
13  GO:1901019 regulation of calcium ion transmembrane ...         7           2     0.09    0.00337 0.28571429
14  GO:0051345 positive regulation of hydrolase activit...        84           5     1.09    0.00475 0.05952381

## p French outgroup
1  GO:0042775 mitochondrial ATP synthesis coupled elec...        36           5     0.60    0.00032 0.13888889
3  GO:0006857                      oligopeptide transport        78           6     1.31    0.00201 0.07692308
4  GO:0009150     purine ribonucleotide metabolic process       462          17     7.76    0.00226 0.03679654
6  GO:0042538              hyperosmotic salinity response        84           6     1.41    0.00292 0.07142857
8  GO:0032787       monocarboxylic acid metabolic process       757          23    12.72    0.00481 0.03038309

#anc French outgroup
1  GO:0046580 negative regulation of Ras protein signa...        15           3     0.11    0.00015 0.20000000
10 GO:0051123 RNA polymerase II preinitiation complex ...        33           3     0.23    0.00160 0.09090909
14 GO:0045944 positive regulation of transcription by ...       156           5     1.10    0.00500 0.03205128

##
go_genes <- genesInTerm(GOdata_xpclr, "GO:0031146")
select_genes[select_genes$V1%in%go_genes$`GO:0031146`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0042538")
select_genes[select_genes$V1%in%go_genes$`GO:0042538`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0018258")
select_genes[select_genes$V1%in%go_genes$`GO:0018258`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0071281")
select_genes[select_genes$V1%in%go_genes$`GO:0071281`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0009414")
s<-select_genes[select_genes$V1%in%go_genes$`GO:0009414`,]
length(unique(s$V1))

allRes_elim[grep(x = allRes_elim$Term,pattern = "water"),]
go_genes <- genesInTerm(GOdata_xpclr, "GO:0009414")
select_genes[select_genes$V1%in%go_genes$`GO:0009414`,]

go_genes <- genesInTerm(GOdata_xpclr, "GO:0009408")
select_genes[select_genes$V1%in%go_genes$`GO:0009408`,]

select_genes[select_genes$V1=="AA3G31360",]
##FR
#PB
5552 AA2G15790 5.161387
9185 AA3G22390 4.227661
5572 AA2G15990 3.663618
#AB
26639 AA8G20000 2.695309 # in AB of W-CAN
2402  AA1G25010 2.295034
15769 AA5G09300 2.485690

#ECAN
#PB
9185   AA3G22390 14.625835
8833   AA3G18870 11.020767, 88331  AA3G18870  5.673879, 88333  AA3G18870  6.416051, 88332  AA3G18870  6.779682
16816  AA5G19770  8.904166
25689  AA8G10500  7.086285
25690  AA8G10510  7.086285
20806  AA6G27810  6.941917,208061 AA6G27810  6.257713
9592   AA3G26460  6.027222
#AB
#none

allRes_elim[grep(x = allRes_elim$Term,pattern = "heat"),]
go_genes <- genesInTerm(GOdata_xpclr, "GO:1900034")
allRes_elim[grep(x = allRes_elim$Term,pattern = "drought"),]
  go_genes <- genesInTerm(GOdata_xpclr, "GO:1900034")

select_genes[select_genes$V1%in%go_genes$`GO:1900034`,]
go_genes <- genesInTerm(GOdata_xpclr, "GO:0010286")
select_genes[select_genes$V1%in%go_genes$`GO:0010286`,]
go_genes <- genesInTerm(GOdata_xpclr, "GO:0009408")
select_genes[select_genes$V1%in%go_genes$`GO:0009408`,]
##FR
#PB
26147 AA8G15080 3.354254
26148 AA8G15090 3.354254
#AB
AA8G20000 2.695309
##ECAN
#PB
20808 AA6G27830 5.37963
#AB
2849 AA1G29480 6.551099


gff_select[order(gff_select$score, decreasing = T),]
gff_select[order(gff_select$score, decreasing = T),]

gff_select[]

allRes_select[ 1:40,]



#
"AA8G47450" %in% gff_select$V9

###in 3 GPAT AA4G30960, RUP1 AA8G30040, VIP4 AA8G47450, EZA1 AA6G21780, VIN3-like AA3G25750
## in 3.1 Anc 1% tail MYB4 AA6G21460, CSP4 AA3G18330, NAC054 AA3G18330, EAF1A AA3G25230
## in 3.2 and 3.3 1% tail EZA1, VIN3-like-1, EAF1A AA3G25230
## in 3.2 A GA20OX2 AA8G29590, VIP4

###in 2 private branches: LSU2 AA8G45080, RUP1 AA8G30040, SFH3 ID=AA6G22060, LOX3 AA1G22960, BAM2 AA5G14930
###in 2 anc GA20OX4 AA2G06360, RUP1 AA8G30040, KNAT2 3.170924

## 2.1 B 1% AA2G15160 APATELA 1










#########plotting
#color set for plotting
set_color<- colors <- c(
  "#1f78b4",  # blue 2.1
  "#e66101",  # orange 2.2
  "#33a02c",  # green 2.3
  "#e41a1c",  # red 3.1
  "#984ea3",  # purple 3.2
  "#ffb000"   # golden yellow 3.3
)
#set shapes for groups
set_shapes <- c(
  21, # circle W-CAN1 (ES03)
  22, # rect W-CAN2 (ES04)
  33 # diamond ES17
)

set_lines <- c(
  1, # circle W-CAN1 (ES03)
  2, # rect W-CAN2 (ES04)
  4 # diamond ES17
)

#NAC055
ggplot(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3",],aes(x=PhysPos, y=X3PCLR.B.std))+
  geom_point()+
  geom_path()+
  scale_x_continuous(limits=c(10030000, 10130000))

ggplot(data = Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr3",],aes(x=PhysPos, y=X3PCLR.B.std))+
  geom_point()+
  geom_path()+
  scale_x_continuous(limits=c(10030000, 10130000))

mask_CN$type.0.single.copy.1.multi.copy <- as.factor(mask_CN$type.0.single.copy.1.multi.copy)
mask_FR$type.0.single.copy.1.multi.copy <- as.factor(mask_FR$type.0.single.copy.1.multi.copy)
mask_CN_4plot <- mask_CN[mask_CN$Chromosome=="chr3" & mask_CN$End>=10030000 & mask_CN$Start <=10130000,]
mask_CN_4plot$Start[1]<- 10030000
mask_CN_4plot$End[nrow(mask_CN_4plot)]<- 10130000
mask_FR_4plot <- mask_FR[mask_FR$Chromosome=="chr3" & mask_FR$End>=10030000 & mask_FR$Start <=10130000,]
mask_FR_4plot$Start[1]<- 10030000
mask_FR_4plot$End[nrow(mask_FR_4plot)]<- 10130000

SetAll_xpclr_chrall$set <- factor(SetAll_xpclr_chrall$set)
# plot_peak_NAC055 <-
ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_line(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3" & Set3.1_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1, colour = set_color[4])+
  geom_line(data = Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr3" & Set3.3_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1, colour = set_color[6])+
  # geom_line(data = Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr3" & Set3.2_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1)+
  # geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std, color=set),size=1.2, alpha=1)+
  # geom_line()+
  scale_y_continuous(limits = c(-2.2,12), expand = c(0.1,0), breaks=seq(-2,12, by=2))+
  scale_x_continuous(limits=c(9970000, 10130000), expand = c(0.1,0), breaks = seq(10030000, 10130000, by=20000), labels=seq(10.03, 10.13, by=0.02))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  #unknown gene
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10047205, xend=10048334,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #trasnposon
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10075217, xend=10076201,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #NAC055
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10077093, xend=10078587,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #NAC056
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10106697, xend=10108574,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+

  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  # scale_color_manual(name="",breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],y=-2.2 - (12+2.2)*0.1, yend=-2.2 - (12+2.2)*0.1, x=10030000, xend=10130000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],y=-2, yend=12, x=10030000+ (10030000 -10130000)*0.1, xend=10030000+ (10030000 -10130000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="Chromosome 3")+
  geom_rect(data = mask_CN_4plot, aes(xmin=Start, xmax=End, ymin=-2.2, ymax=-1.9, fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  geom_rect(data = mask_FR_4plot, aes(xmin=Start, xmax=End, ymin=11.7, ymax=12,fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=40), axis.title=element_text(size=40), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=40, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  annotate(geom = "text",x =  10047205+(10048334 -10047205)/2, y = -0.3, label=expression(italic("AA3G18860")), size=8, parse=T)+
  annotate(geom = "text",x =  10075217+(10076201 -10075217)/2, y = -1.5, label=expression(italic("hAT-Ac")), size=8, parse=T)+
  annotate(geom = "text",x =  10077093+(10078587 -10077093)/2, y = -0.3, label=expression(italic("NAC055")), size=8, parse=T)+
  annotate(geom = "text",x =  10106697+(10108574 -10106697)/2, y = -0.3, label=expression(italic("NAC056")), size=8, parse=T)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))



mask_CN$type.0.single.copy.1.multi.copy <- as.factor(mask_CN$type.0.single.copy.1.multi.copy)
mask_FR$type.0.single.copy.1.multi.copy <- as.factor(mask_FR$type.0.single.copy.1.multi.copy)
mask_CN_4plot <- mask_CN[mask_CN$Chromosome=="chr3" & mask_CN$End>=9970000 & mask_CN$Start <=10130000,]
mask_CN_4plot$Start[1]<- 9970000
mask_CN_4plot$End[nrow(mask_CN_4plot)]<- 10130000
mask_FR_4plot <- mask_FR[mask_FR$Chromosome=="chr3" & mask_FR$End>=9970000 & mask_FR$Start <=10130000,]
mask_FR_4plot$Start[1]<- 9970000
mask_FR_4plot$End[nrow(mask_FR_4plot)]<- 10130000

SetAll_xpclr_chrall$set <- factor(SetAll_xpclr_chrall$set)
# plot_peak_NAC055 <-
p1<- ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_line(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3" & Set3.1_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1, colour = set_color[4])+
  geom_line(data = Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr3" & Set3.3_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1, colour = set_color[6])+
  # geom_line(data = Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr3" & Set3.2_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1)+
  # geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std, color=set),size=1.2, alpha=1)+
  # geom_line()+
  scale_y_continuous(limits = c(-2.2,14), expand = c(0.1,0), breaks=seq(-2,14, by=2))+
  scale_x_continuous(limits=c(9970000, 10130000), expand = c(0.1,0), breaks = seq(9970000, 10130000, by=40000), labels=seq(9.97, 10.13, by=0.04))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  #unknown gene
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=9971744, xend=9975491,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=9993856, xend=9995638,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10008126, xend=10011788,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=1.4, yend=1.4, x=10012333, xend=10013691,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=0.2, yend=0.2, x=10014619, xend=10016741,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10018411, xend=10017244,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10022010, xend=10022120,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #unknown gene
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10047205, xend=10048334,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #trasnposon
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10075217, xend=10076201,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #NAC055
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10077093, xend=10078587,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #NAC056
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10106697, xend=10108574,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+

  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  # scale_color_manual(name="",breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],y=-2.2 - (14+2.2)*0.1, yend=-2.2 - (14+2.2)*0.1, x=9970000, xend=10130000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],y=-2, yend=14, x=9970000+ (9970000 -10130000)*0.1, xend=9970000+ (9970000 -10130000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  geom_rect(data = mask_CN_4plot, aes(xmin=Start, xmax=End, ymin=-2.2, ymax=-1.9, fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  geom_rect(data = mask_FR_4plot, aes(xmin=Start, xmax=End, ymin=13.7, ymax=14,fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=40), axis.title=element_text(size=40), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=40, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  annotate(geom = "text",x =  9971744+(9975491 -9971744)/2, y = -0.3, label=expression(italic("CKL8")), size=6, parse=T)+
  annotate(geom = "text",x =  9993856+(9995638 -9993856)/2, y = -0.3, label=expression(italic("AA3G18800")), size=6, parse=T)+
  annotate(geom = "text",x =  10008126+(10011788 -10008126)/2, y = -0.3, label=expression(italic("WDR44")), size=6, parse=T)+
  annotate(geom = "text",x =  10012333+(10013691 -10012333)/2, y =  2.1, label=expression(italic("AA3G18820")), size=6, parse=T)+
  annotate(geom = "text",x =  10014619+(10016741 -10014619)/2, y =  0.9, label=expression(italic("BRIX1-1")), size=6, parse=T)+
  annotate(geom = "text",x =  10017244+(10018411 -10017244)/2, y = -1.5, label=expression(italic("TSJT1")), size=6, parse=T)+
  annotate(geom = "text",x =  10022010+(10022120 -10022010)/2, y = -0.3, label=expression(italic("AA3G18850")), size=6, parse=T)+
  annotate(geom = "text",x =  10047205+(10048334 -10047205)/2, y = -0.3, label=expression(italic("AA3G18860")), size=6, parse=T)+
  annotate(geom = "text",x =  10075217+(10076201 -10075217)/2, y = -1.5, label=expression(italic("hAT-Ac")), size=6, parse=T)+
  annotate(geom = "text",x =  10077093+(10078587 -10077093)/2, y = -0.3, label=expression(italic("NAC055")), size=6, parse=T)+
  annotate(geom = "text",x =  10106697+(10108574 -10106697)/2, y = -0.3, label=expression(italic("NAC056")), size=6, parse=T)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))


p2<-ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3" ,],aes(x=PhysPos, y=X3PCLR.B.std),size=3, alpha=0.5, shape=21, color="black",fill = set_color[4])+
  geom_point(data = Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr3" ,],aes(x=PhysPos, y=X3PCLR.B.std),size=3, alpha=0.5, shape=21, color="black", fill = set_color[6])+
  # geom_line()+
  scale_y_continuous(limits = c(-2.2,14), expand = c(0.1,0), breaks=seq(-2,14, by=2))+
  scale_x_continuous(limits=c(0, 30421000), expand = c(0.1,0),breaks = seq(0, 30000000, by=10e6), labels=seq(0, 30, by=10))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (14+2.2)*0.1, yend=-2.2 - (14+2.2)*0.1, x=0, xend=30000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=14, x=-30421000*0.1, xend=-30421000*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=40), axis.title=element_text(size=40), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=40), legend.title = element_text(size=40), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=40, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))


ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$set %in% c("3.1", "3.3") ,],aes(x=PhysPos, y=X3PCLR.B.std, fill=set),size=3, alpha=1, shape=21, color="black")+
  geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$set %in% c("3.1", "3.3") ,],aes(x=PhysPos, y=X3PCLR.B.std, color =  set),size=1.2)+
  # geom_line()+
  scale_y_continuous(limits = c(-2.2,14), expand = c(0.1,0), breaks=seq(-2,14, by=2))+
  scale_x_continuous(limits=c(0, 30421000), expand = c(0.1,0),breaks = seq(0, 30000000, by=10e6), labels=seq(0, 30, by=10))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  scale_color_manual(breaks = c("3.1","3.3"), values = set_color[c(4,6)], labels=c("Set4: P(W-CAN1)","Set6: P(ES17)"))+
  scale_fill_manual(breaks = c("3.1","3.3"), values = set_color[c(4,6)],labels=c("Set4: P(W-CAN1)","Set6: P(ES17)"))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (14+2.2)*0.1, yend=-2.2 - (14+2.2)*0.1, x=0, xend=30000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=14, x=-30421000*0.1, xend=-30421000*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="", color="", fill="")+
  theme(panel.background = element_blank(), axis.text = element_text(size=40), axis.title=element_text(size=40), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=20), legend.title = element_text(size=10), legend.key.size = unit(1,"cm") , legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=40, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, size=5)))


library(cowplot)
plot_grid(plotlist = list(p1,p2), align = "hv", nrow = 2)






#### flowering genes
##VIP4 chr8    maker   gene    46056760        46060745

SetAll_xpclr_chrall$set <- as.factor(SetAll_xpclr_chrall$set)
ggplot(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],aes(x=PhysPos, y=X3PCLR.A.std, color=set))+
  geom_point()+
  geom_path()+
  scale_x_continuous(limits=c(46000000, 46120000))


mask_CN$type.0.single.copy.1.multi.copy <- as.factor(mask_CN$type.0.single.copy.1.multi.copy)
mask_FR$type.0.single.copy.1.multi.copy <- as.factor(mask_FR$type.0.single.copy.1.multi.copy)
mask_CN_4plot <- mask_CN[mask_CN$Chromosome=="chr8" & mask_CN$End>=46000000 & mask_CN$Start <=46120000,]
mask_CN_4plot$Start[1]<- 46000000
mask_CN_4plot$End[nrow(mask_CN_4plot)]<- 46120000
mask_FR_4plot <- mask_FR[mask_FR$Chromosome=="chr8" & mask_FR$End>=46000000 & mask_FR$Start <=46120000,]
mask_FR_4plot$Start[1]<- 46000000
mask_FR_4plot$End[nrow(mask_FR_4plot)]<- 46120000




SetAll_xpclr_chrall$set <- factor(SetAll_xpclr_chrall$set)
# plot_peak_NAC055 <-

c3.1A<-sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]

p_vip4 <- ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_line(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr8" & Set3.1_xpclr_chrall$PhysPos<=46120000,],aes(x=PhysPos, y=X3PCLR.A.std),size=2, alpha=1, colour = set_color[4])+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[4], y=c3.1A, yend=c3.1A, x=46000000, xend=46120000, linetype = 3, size=2)+
  # geom_line(data = Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr3" & Set3.2_xpclr_chrall$PhysPos<=46120000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1)+
  # geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$PhysPos<=46120000,],aes(x=PhysPos, y=X3PCLR.B.std, color=set),size=1.2, alpha=1)+
  # geom_line()+
  scale_y_continuous(limits = c(-4.4,24), expand = c(0.1,0), breaks=seq(-4,24, by=4))+
  scale_x_continuous(limits=c(46000000, 46120000), expand = c(0.1,0), breaks = seq(46000000, 46120000, by=60000), labels=seq(46, 46.12, by=0.06))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  #unknown gene
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=46042942, xend=46040146,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=1.5, yend=1.5, x=46043497, xend=46056637,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=46060745, xend=46056760,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  # scale_color_manual(name="",breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],y=-4.4 - (24+4.4)*0.1, yend=-4.4 - (24+4.4)*0.1, x=46000000, xend=46120000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],y=-4, yend=24, x=46000000+ (46000000 -46120000)*0.1, xend=46000000+ (46000000 -46120000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  geom_rect(data = mask_CN_4plot, aes(xmin=Start, xmax=End, ymin=-4.4, ymax=-3.9, fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  geom_rect(data = mask_FR_4plot, aes(xmin=Start, xmax=End, ymin=23.5, ymax=24,fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=45), axis.title=element_text(size=45), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=45), legend.title = element_text(size=45), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=45, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  annotate(geom = "text",x =  46040146+(46042942 -46040146)/2, y = 0, label=expression(italic("PDCB1")), size=8, parse=T)+
  annotate(geom = "text",x =  46043497+(46056637 -46043497)/2, y = 2.5, label=expression(italic("BRRC2")), size=8, parse=T)+
  annotate(geom = "text",x =  46056760+(46060745 -46056760)/2, y = 0, label=expression(italic("VIP4")), size=8, parse=T)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))


# p2<-
p_chr8_S3.1<-ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr8" ,],aes(x=PhysPos, y=X3PCLR.A.std),size=5, alpha=0.5, shape=21, color="black",fill = set_color[4])+
  # geom_line()+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[4], y=c3.1A, yend=c3.1A, x=0, xend=50000000, linetype = 3, size=2)+
  scale_y_continuous(limits = c(-4.4,24), expand = c(0.1,0), breaks=seq(-4,24, by=4))+
  scale_x_continuous(limits=c(0, 50000000), expand = c(0.1,0),breaks = seq(0, 50000000, by=10e6), labels=seq(0, 50, by=10))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],y=-4.4 - (24+4.4)*0.1, yend=-4.4 - (24+4.4)*0.1, x= 0, xend=50000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],y=-4, yend=24, x=0+ (0 -50000000)*0.1, xend=0+ (0 -50000000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=45), axis.title=element_text(size=45), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=45), legend.title = element_text(size=45), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=45, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))


tail(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr8",])
max(Set3.1_xpclr_chrall$X3PCLR.A.std)


ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$set %in% c("3.1") ,],aes(x=PhysPos, y=X3PCLR.B.std, fill=set),size=3, alpha=1, shape=21, color="black")+
  geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$set %in% c("3.1") ,],aes(x=PhysPos, y=X3PCLR.B.std, color =  set),size=1.2)+
  # geom_line()+
  scale_y_continuous(limits = c(-2.2,14), expand = c(0.1,0), breaks=seq(-2,14, by=2))+
  scale_x_continuous(limits=c(0, 30421000), expand = c(0.1,0),breaks = seq(0, 30000000, by=10e6), labels=seq(0, 30, by=10))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  scale_color_manual(breaks = c("3.1"), values = set_color[c(4)], labels=c("Set4: P(W-CAN2)"))+
  scale_fill_manual(breaks = c("3.1"), values = set_color[c(4)],labels=c("Set4: P(W-CAN2)"))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (14+2.2)*0.1, yend=-2.2 - (14+2.2)*0.1, x=0, xend=30000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=14, x=-30421000*0.1, xend=-30421000*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="", color="", fill="")+
  theme(panel.background = element_blank(), axis.text = element_text(size=40), axis.title=element_text(size=40), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=20), legend.title = element_text(size=10), legend.key.size = unit(1,"cm") , legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=40, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, size=5)))





##VIL1 chr3    maker   gene    17794108        17796816

SetAll_xpclr_chrall$set <- as.factor(SetAll_xpclr_chrall$set)
ggplot(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],aes(x=PhysPos, y=X3PCLR.B.std, color=set))+
  geom_point()+
  geom_path()

ggplot(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],aes(x=PhysPos, y=X3PCLR.B.std, color=set))+
  geom_point()+
  geom_path()+
  scale_x_continuous(limits=c(17650000, 17830000))

ggplot(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],aes(x=PhysPos, y=X3PCLR.B.std, color=set))+
  geom_point()+
  geom_path()+
  scale_x_continuous(limits=c(17650000, 17830000))



mask_CN$type.0.single.copy.1.multi.copy <- as.factor(mask_CN$type.0.single.copy.1.multi.copy)
mask_FR$type.0.single.copy.1.multi.copy <- as.factor(mask_FR$type.0.single.copy.1.multi.copy)
mask_CN_4plot <- mask_CN[mask_CN$Chromosome=="chr3" & mask_CN$End>=17650000 & mask_CN$Start <=17830000,]
mask_CN_4plot$Start[1]<- 17650000
mask_CN_4plot$End[nrow(mask_CN_4plot)]<- 17830000
mask_FR_4plot <- mask_FR[mask_FR$Chromosome=="chr3" & mask_FR$End>=17650000 & mask_FR$Start <=17830000,]
mask_FR_4plot$Start[1]<- 17650000
mask_FR_4plot$End[nrow(mask_FR_4plot)]<- 17830000

SetAll_xpclr_chrall$set <- factor(SetAll_xpclr_chrall$set)
# plot_peak_NAC055 <-

c3.2B<-sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.B.std)]
c3.3B<-sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set3.3_xpclr_chrall$X3PCLR.B.std)]

p_VIL1_S3.3<-ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  # geom_line(data = Set2.2_xpclr_chrall[Set2.2_xpclr_chrall$Chr=="chr3" & Set2.2_xpclr_chrall$PhysPos<=17830000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1, colour = set_color[2])+
  # geom_line(data = Set2.3_xpclr_chrall[Set2.3_xpclr_chrall$Chr=="chr3" & Set2.3_xpclr_chrall$PhysPos<=17830000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1, colour = set_color[3])+
  # geom_line(data = Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr3" & Set3.2_xpclr_chrall$PhysPos<=17830000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1, colour = set_color[5])+
  geom_line(data = Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr3" & Set3.3_xpclr_chrall$PhysPos<=17830000,],aes(x=PhysPos, y=X3PCLR.B.std),size=2, alpha=1, colour = set_color[6])+
  # geom_hline(yintercept = c3.3B, color=set_color[6], linetype = "dashed")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[6], y=c3.3B, yend=c3.3B, x=17650000, xend=17830000, linetype = 3, size=2)+
  # geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[5], y=c3.2B-0.1, yend=c3.2B, x=17650000, xend=17830000, linetype = 2)+
  # geom_line(data = Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr3" & Set3.2_xpclr_chrall$PhysPos<=17830000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1)+
  # geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$PhysPos<=17830000,],aes(x=PhysPos, y=X3PCLR.B.std, color=set),size=1.2, alpha=1)+
  # geom_line()+
  scale_y_continuous(limits = c(-2.2,12), expand = c(0.1,0), breaks=seq(-2,12, by=2))+
  scale_x_continuous(limits=c(17650000, 17830000), expand = c(0.1,0), breaks = seq(17650000, 17830000, by=90000), labels=seq(17.65, 17.83, by=0.09))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  #unknown gene
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=17681360, xend=17686475,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=0.5, yend=0.5, x=17690830, xend=17686542,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-0.25, yend=-0.25, x=17770301, xend=17772349,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=0.5, yend=0.5, x=17794108, xend=17796816,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=17799583, xend=17798414,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  # scale_color_manual(name="",breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (12+2.2)*0.1, yend=-2.2 - (12+2.2)*0.1, x=17650000, xend=17830000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=12, x=17650000+ (17650000 -17830000)*0.1, xend=17650000+ (17650000 -17830000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  geom_rect(data = mask_CN_4plot, aes(xmin=Start, xmax=End, ymin=-2.2, ymax=-1.9, fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  geom_rect(data = mask_FR_4plot, aes(xmin=Start, xmax=End, ymin=11.7, ymax=12,fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=45), axis.title=element_text(size=45), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=45), legend.title = element_text(size=45), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=45, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  annotate(geom = "text",x =  17681360+(17686475 -17681360)/2, y = -0.3, label=expression(italic("AA3G25720")), size=8, parse=T)+
  annotate(geom = "text",x =  17690830+(17686542 -17690830)/2, y = 1.2, label=expression(italic("AA3G25730")), size=8, parse=T)+
  annotate(geom = "text",x =  17770301+(17772349 -17770301)/2, y = 0.45, label=expression(italic("NAKR1")), size=8, parse=T)+
  annotate(geom = "text",x =  17794108+(17796816 -17794108)/2, y = 1.2, label=expression(italic("VIL1")), size=8, parse=T)+
  annotate(geom = "text",x =  17798414+(17799583 -17798414)/2, y = -0.3, label=expression(italic("AA3G25760")), size=8, parse=T)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))

# chr3    maker   gene    17794108        17796816
# chr3    maker   gene    17798414        17799583 AA3G25760
  # p2<-
p_chr3_S3.3<- ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr3" ,],aes(x=PhysPos, y=X3PCLR.B.std),size=5, alpha=0.5, shape=21, color="black",fill = set_color[6])+
  # geom_line()+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[6], y=c3.3B, yend=c3.3B, x=0, xend=30421000, linetype = 3, size=2)+
  scale_y_continuous(limits = c(-2.2,12), expand = c(0.1,0), breaks=seq(-2,12, by=2))+
  scale_x_continuous(limits=c(0, 30421000), expand = c(0.1,0),breaks = seq(0, 50000000, by=10e6), labels=seq(0, 50, by=10))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (12+2.2)*0.1, yend=-2.2 - (12+2.2)*0.1, x= 0, xend=30000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=12, x=0+ (0 -30421000)*0.1, xend=0+ (0 -30421000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=45), axis.title=element_text(size=45), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=45), legend.title = element_text(size=45), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=45, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))


tail(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3",])
max(Set3.1_xpclr_chrall$X3PCLR.A.std)


ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$set %in% c("3.3") ,],aes(x=PhysPos, y=X3PCLR.B.std, fill=set),size=3, alpha=1, shape=21, color="black")+
  geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$set %in% c("3.3") ,],aes(x=PhysPos, y=X3PCLR.B.std, color =  set),size=1.2)+
  # geom_line()+
  scale_y_continuous(limits = c(-2.2,14), expand = c(0.1,0), breaks=seq(-2,14, by=2))+
  scale_x_continuous(limits=c(0, 30421000), expand = c(0.1,0),breaks = seq(0, 30000000, by=10e6), labels=seq(0, 30, by=10))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  scale_color_manual(breaks = c("3.3"), values = set_color[c(6)], labels=c("Set6: P(ES17)"))+
  scale_fill_manual(breaks = c("3.3"), values = set_color[c(6)],labels=c("Set6: P(ES17)"))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (14+2.2)*0.1, yend=-2.2 - (14+2.2)*0.1, x=0, xend=30000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=14, x=-30421000*0.1, xend=-30421000*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="", color="", fill="")+
  theme(panel.background = element_blank(), axis.text = element_text(size=40), axis.title=element_text(size=40), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=20), legend.title = element_text(size=10), legend.key.size = unit(1,"cm") , legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=40, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, size=5)))

##FIP1 chr3    maker   gene    24873902        24878878

c2.1A<-sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]

SetAll_xpclr_chrall$set <- as.factor(SetAll_xpclr_chrall$set)
ggplot(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],aes(x=PhysPos, y=X3PCLR.A.std, color=set))+
  geom_point()+
  geom_path()

ggplot(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],aes(x=PhysPos, y=X3PCLR.A.std, color=set))+
  geom_point()+
  geom_path()+
  scale_x_continuous(limits=c(24850000, 24900000))




mask_CN$type.0.single.copy.1.multi.copy <- as.factor(mask_CN$type.0.single.copy.1.multi.copy)
mask_FR$type.0.single.copy.1.multi.copy <- as.factor(mask_FR$type.0.single.copy.1.multi.copy)
mask_CN_4plot <- mask_CN[mask_CN$Chromosome=="chr3" & mask_CN$End>=24850000 & mask_CN$Start <=24900000,]
mask_CN_4plot$Start[1]<- 24850000
mask_CN_4plot$End[nrow(mask_CN_4plot)]<- 24900000
mask_FR_4plot <- mask_FR[mask_FR$Chromosome=="chr3" & mask_FR$End>=24850000 & mask_FR$Start <=24900000,]
mask_FR_4plot$Start[1]<- 24850000
mask_FR_4plot$End[nrow(mask_FR_4plot)]<- 24900000

SetAll_xpclr_chrall$set <- factor(SetAll_xpclr_chrall$set)
# plot_peak_NAC055 <-
p_FIP1_S3.1_S2.1<-ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_line(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3" & Set3.1_xpclr_chrall$PhysPos<=24900000,],aes(x=PhysPos, y=X3PCLR.A.std),size=2, alpha=1, colour = set_color[4])+
  geom_line(data = Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr3" & Set2.1_xpclr_chrall$PhysPos<=24900000,],aes(x=PhysPos, y=X3PCLR.A.std),size=2, alpha=1, colour = set_color[1])+
  # geom_line(data = Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr3" & Set3.2_xpclr_chrall$PhysPos<=24900000,],aes(x=PhysPos, y=X3PCLR.A.std),size=1.2, alpha=1)+
  # geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$PhysPos<=24900000,],aes(x=PhysPos, y=X3PCLR.A.std, color=set),size=1.2, alpha=1)+
  # geom_line()+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[4], y=c3.1A, yend=c3.1A, x=24850000, xend=24900000, linetype = 3, size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[1], y=c2.1A, yend=c2.1A, x=24850000, xend=24900000, linetype = 3, size=2)+
  scale_y_continuous(limits = c(-2.2,12), expand = c(0.1,0), breaks=seq(-2,12, by=2))+
  scale_x_continuous(limits=c(24850000, 24900000), expand = c(0.1,0), breaks = seq(24850000, 24900000, by=25000), labels=seq(24.85, 24.9, by=0.025))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  #unknown gene
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-0.25, yend=-0.25, x=24871486, xend=24869888,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=24873902, xend=24878878,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=0.5, yend=0.5, x=24882120, xend=24883737,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=24887506, xend=24884256,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  # scale_color_manual(name="",breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (12+2.2)*0.1, yend=-2.2 - (12+2.2)*0.1, x=24850000, xend=24900000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=12, x=24850000+ (24850000 -24900000)*0.1, xend=24850000+ (24850000 -24900000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  geom_rect(data = mask_CN_4plot, aes(xmin=Start, xmax=End, ymin=-2.2, ymax=-1.9, fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  geom_rect(data = mask_FR_4plot, aes(xmin=Start, xmax=End, ymin=11.7, ymax=12,fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=45), axis.title=element_text(size=45), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=45), legend.title = element_text(size=45), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=45, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  annotate(geom = "text",x =  24869888+(24871486 -24869888)/2, y = 0.45, label=expression(italic("SOAR1")), size=8, parse=T)+
  annotate(geom = "text",x =  24873902+(24878878 -24873902)/2, y = -0.3, label=expression(italic("FIP1")), size=8, parse=T)+
  annotate(geom = "text",x =  24882120+(24883737 -24882120)/2, y = 1.2, label=expression(italic("AA3G31370")), size=8, parse=T)+
  annotate(geom = "text",x =  24884256+(24887506 -24884256)/2, y = -0.3, label=expression(italic("FBLX13")), size=8, parse=T)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))

# chr3    maker   gene    17794108        17796816
# chr3    maker   gene    17798414        17799583 AA3G25760
# p2<-
p_chr3_S3.1_S2.1<-ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3" ,],aes(x=PhysPos, y=X3PCLR.A.std),size=5, alpha=0.5, shape=21, color="black",fill = set_color[4])+
  geom_point(data = Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr3" ,],aes(x=PhysPos, y=X3PCLR.A.std),size=5, alpha=0.5, shape=21, color="black",fill = set_color[1])+
  # geom_line()+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[4], y=c3.1A, yend=c3.1A, x=0, xend=30421000, linetype = 3, size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[1], y=c2.1A, yend=c2.1A, x=0, xend=30421000, linetype = 3, size=2)+
  scale_y_continuous(limits = c(-2.2,12), expand = c(0.1,0), breaks=seq(-2,12, by=2))+
  scale_x_continuous(limits=c(0, 30421000), expand = c(0.1,0),breaks = seq(0, 50000000, by=10e6), labels=seq(0, 50, by=10))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (12+2.2)*0.1, yend=-2.2 - (12+2.2)*0.1, x= 0, xend=30000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=12, x=0+ (0 -30421000)*0.1, xend=0+ (0 -30421000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=45), axis.title=element_text(size=45), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=45), legend.title = element_text(size=45), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=45, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))

tail(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3",])
max(Set3.1_xpclr_chrall$X3PCLR.A.std)
min(Set3.1_xpclr_chrall$X3PCLR.A.std[Set3.1_xpclr_chrall$Chr=="chr3"])
min(Set2.1_xpclr_chrall$X3PCLR.A.std[Set2.1_xpclr_chrall$Chr=="chr3"])







library(cowplot)
p_flowering<-plot_grid(plotlist = list(p_vip4, p_VIL1_S3.3, p_FIP1_S3.1_S2.1, p_chr8_S3.1, p_chr3_S3.3, p_chr3_S3.1_S2.1), align = "hv", nrow = 2)
p_flowering

ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$set %in% c("2.1","3.1") ,],aes(x=PhysPos, y=X3PCLR.B.std, fill=set),size=3, alpha=1, shape=21, color="black")+
  geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$set %in% c("2.1","3.1") ,],aes(x=PhysPos, y=X3PCLR.B.std, color =  set),size=1.2)+
  # geom_line()+
  scale_y_continuous(limits = c(-2.2,14), expand = c(0.1,0), breaks=seq(-2,14, by=2))+
  scale_x_continuous(limits=c(0, 30421000), expand = c(0.1,0),breaks = seq(0, 30000000, by=10e6), labels=seq(0, 30, by=10))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  scale_color_manual(breaks = c("2.1","3.1"), values = set_color[c(1,4)], labels=c("Set1: P(W-CAN2)","Set4: P(W-CAN2)"))+
  scale_fill_manual(breaks = c("2.1","3.1"), values = set_color[c(1,4)],labels=c("Set1: P(W-CAN2)","Set4: P(W-CAN2)"))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (14+2.2)*0.1, yend=-2.2 - (14+2.2)*0.1, x=0, xend=30000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=14, x=-30421000*0.1, xend=-30421000*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="", color="", fill="")+
  theme(panel.background = element_blank(), axis.text = element_text(size=40), axis.title=element_text(size=40), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=20), legend.title = element_text(size=10), legend.key.size = unit(1,"cm") , legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=40, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, size=5)))


###other high peaks SC5D

## AA6G27580 chr6    maker   gene    28457738        28459744

c2.1A<-sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]
c3.1A<-sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]
c3.2A<-sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]

c3.1B<-sort(Set3.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.B.std)]
c3.3B<-sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set3.3_xpclr_chrall$X3PCLR.B.std)]



SetAll_xpclr_chrall$set <- as.factor(SetAll_xpclr_chrall$set)
ggplot(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr6",],aes(x=PhysPos, y=X3PCLR.A.std, color=set))+
  geom_point()+
  geom_path()

ggplot(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr6",],aes(x=PhysPos, y=X3PCLR.A.std, color=set))+
  geom_point()+
  geom_path()+
  scale_x_continuous(limits=c(28250000, 28550000))




mask_CN$type.0.single.copy.1.multi.copy <- as.factor(mask_CN$type.0.single.copy.1.multi.copy)
mask_FR$type.0.single.copy.1.multi.copy <- as.factor(mask_FR$type.0.single.copy.1.multi.copy)
mask_CN_4plot <- mask_CN[mask_CN$Chromosome=="chr6" & mask_CN$End>=28250000 & mask_CN$Start <=28550000,]
mask_CN_4plot$Start[1]<- 28250000
mask_CN_4plot$End[nrow(mask_CN_4plot)]<- 28550000
mask_FR_4plot <- mask_FR[mask_FR$Chromosome=="chr6" & mask_FR$End>=28250000 & mask_FR$Start <=28550000,]
mask_FR_4plot$Start[1]<- 28250000
mask_FR_4plot$End[nrow(mask_FR_4plot)]<- 28550000

SetAll_xpclr_chrall$set <- factor(SetAll_xpclr_chrall$set)
# plot_peak_NAC055 <-
p_SC5D<-ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_line(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr6" & Set3.1_xpclr_chrall$PhysPos<=28550000,],aes(x=PhysPos, y=X3PCLR.A.std),size=2, alpha=1, colour = set_color[4])+
  geom_line(data = Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr6" & Set3.2_xpclr_chrall$PhysPos<=28550000,],aes(x=PhysPos, y=X3PCLR.A.std),size=2, alpha=1, colour = set_color[5])+
  geom_line(data = Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr6" & Set2.1_xpclr_chrall$PhysPos<=28550000,],aes(x=PhysPos, y=X3PCLR.A.std),size=2, alpha=1, colour = set_color[1])+
  # geom_line(data = Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr6" & Set3.2_xpclr_chrall$PhysPos<=28550000,],aes(x=PhysPos, y=X3PCLR.A.std),size=1.2, alpha=1)+
  # geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr6" & SetAll_xpclr_chrall$PhysPos<=28550000,],aes(x=PhysPos, y=X3PCLR.A.std, color=set),size=1.2, alpha=1)+
  # geom_line()+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[4], y=c3.1A, yend=c3.1A, x=28250000, xend=28550000, linetype = 3, size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[5], y=c3.2A, yend=c3.2A, x=28250000, xend=28550000, linetype = 3, size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[1], y=c2.1A, yend=c2.1A, x=28250000, xend=28550000, linetype = 3, size=2)+
  scale_y_continuous(limits = c(-4.4,20), expand = c(0.1,0), breaks=seq(-4,20, by=4))+
  scale_x_continuous(limits=c(28250000, 28550000), expand = c(0.1,0), breaks = seq(28250000, 28550000, by=100000), labels=seq(28.25, 28.55, by=0.1))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  #unknown gene
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-2.5, yend=-2.5, x=28272269, xend=28274927,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=8, yend=8, x=28275827, xend=28275243,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y= 12, yend=12, x=28278221, xend=28276445,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=10, yend=10, x=28286344, xend=28279302,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1.5, yend=-1.5, x=28310995, xend=28309191,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=28375166, xend=28373717,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-2.5, yend=-2.5, x=28450479, xend=28447626,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-0.5, yend=-0.5, x=28459744, xend=28457738,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  # scale_color_manual(name="",breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr6",],y=-4.4 - (20+4.4)*0.1, yend=-4.4 - (20+4.4)*0.1, x=28250000, xend=28550000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr6",],y=-4, yend=20, x=28250000+ (28250000 -28550000)*0.1, xend=28250000+ (28250000 -28550000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  geom_rect(data = mask_CN_4plot, aes(xmin=Start, xmax=End, ymin=-4.4, ymax=-3.6, fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  geom_rect(data = mask_FR_4plot, aes(xmin=Start, xmax=End, ymin=19.2, ymax=20,fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=45), axis.title=element_text(size=45), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=45), legend.title = element_text(size=40), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,0,1), "cm"),  plot.title = element_text(size=45, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  annotate(geom = "text",x =  28272269+(28274927 -28272269)/2, y = -1.5, label=expression(italic("STP7")), size=8, parse=T)+
  annotate(geom = "text",x =  28275243+(28275827 -28275243)/2, y = 9, label=expression(italic("AA6G27520")), size=8, parse=T)+
  annotate(geom = "text",x =  28276445+(28278221 -28276445)/2, y = 13, label=expression(italic("PCMP-E37")), size=8, parse=T)+
  annotate(geom = "text",x =  28279302+(28286344 -28279302)/2, y = 11, label=expression(italic("AA6G27540")), size=8, parse=T)+
  annotate(geom = "text",x =  28309191+(28310995 -28309191)/2, y = -2.5, label=expression(italic("AA6G27550")), size=8, parse=T)+
  annotate(geom = "text",x =  28373717+(28375166 -28373717)/2, y = 0, label=expression(italic("AA6G27560")), size=8, parse=T)+
  annotate(geom = "text",x =  28447626+(28450479 -28447626)/2, y = -1.5, label=expression(italic("GLN1-1")), size=8, parse=T)+
  annotate(geom = "text",x =  28457738+(28459744 -28457738)/2, y = 0.5, label=expression(italic("SC5D")), size=8, parse=T)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))

# chr6    maker   gene    28257823        28260598
# chr6    maker   gene    28275243        28275827
# chr6    maker   gene    28276445        28278221
# chr6    maker   gene    28279302        28286344
# p2<-
p_chr6_S3.1_S2.1<-ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr6" ,],aes(x=PhysPos, y=X3PCLR.A.std),size=5, alpha=0.5, shape=21, color="black",fill = set_color[4])+
  geom_point(data = Set2.1_xpclr_chrall[Set2.1_xpclr_chrall$Chr=="chr6" ,],aes(x=PhysPos, y=X3PCLR.A.std),size=5, alpha=0.5, shape=21, color="black",fill = set_color[1])+
  geom_point(data = Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr6" ,],aes(x=PhysPos, y=X3PCLR.A.std),size=5, alpha=0.5, shape=21, color="black",fill = set_color[5])+
  # geom_line()+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[4], y=c3.1A, yend=c3.1A, x=0, xend=33360400, linetype = 3, size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[5], y=c3.2A, yend=c3.2A, x=0, xend=33360400, linetype = 3, size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[1], y=c2.1A, yend=c2.1A, x=0, xend=33360400, linetype = 3, size=2)+
  scale_y_continuous(limits = c(-4.4,20), expand = c(0.1,0), breaks=seq(-4,20, by=4))+
  scale_x_continuous(limits=c(0, 35000000), expand = c(0.1,0),breaks = seq(0, 35000000, by=5e6), labels=seq(0, 35, by=5))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr6",],y=-4.4 - (20+4.4)*0.1, yend=-4.4 - (20+4.4)*0.1, x= 0, xend=35000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr6",],y=-4, yend=20, x=0+ (0 -35000000)*0.1, xend=0+ (0 -35000000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=45), axis.title=element_text(size=45), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=45), legend.title = element_text(size=45), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(0,1,1,1), "cm"),  plot.title = element_text(size=45, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))

tail(Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr6",])
max(Set3.1_xpclr_chrall$X3PCLR.A.std)
min(Set3.1_xpclr_chrall$X3PCLR.A.std[Set3.1_xpclr_chrall$Chr=="chr6"])
min(Set2.1_xpclr_chrall$X3PCLR.A.std[Set2.1_xpclr_chrall$Chr=="chr6"])


# library(cowplot)
# plot_grid(plotlist = list(p_SC5D, p_chr6_S3.1_S2.1), nrow = 2, align = "hv")


ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$set %in% c("2.1","3.1", "3.2") ,],aes(x=PhysPos, y=X3PCLR.B.std, fill=set),size=3, alpha=1, shape=21, color="black")+
  geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$set %in% c("2.1","3.1", "3.2") ,],aes(x=PhysPos, y=X3PCLR.B.std, color =  set),size=1.2)+
  # geom_line()+
  scale_y_continuous(limits = c(-2.2,14), expand = c(0.1,0), breaks=seq(-2,14, by=2))+
  scale_x_continuous(limits=c(0, 30421000), expand = c(0.1,0),breaks = seq(0, 30000000, by=10e6), labels=seq(0, 30, by=10))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  scale_color_manual(breaks = c("2.1","3.1","3.2"), values = set_color[c(1,4,5)], labels=c("Set1: P(W-CAN2)", "Set4: P(W-CAN2)","Set5: P(W-CAN2)"))+
  scale_fill_manual(breaks = c("2.1","3.1","3.2"), values = set_color[c(1,4,5)],labels=c("Set1: P(W-CAN2)", "Set4: P(W-CAN2)","Set5: P(W-CAN2)"))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (14+2.2)*0.1, yend=-2.2 - (14+2.2)*0.1, x=0, xend=30000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=14, x=-30421000*0.1, xend=-30421000*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="", color="", fill="")+
  theme(panel.background = element_blank(), axis.text = element_text(size=40), axis.title=element_text(size=40), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=20), legend.title = element_text(size=10), legend.key.size = unit(1,"cm") , legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"),  plot.title = element_text(size=40, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, size=5)))

### add NAC plot?
mask_CN$type.0.single.copy.1.multi.copy <- as.factor(mask_CN$type.0.single.copy.1.multi.copy)
mask_FR$type.0.single.copy.1.multi.copy <- as.factor(mask_FR$type.0.single.copy.1.multi.copy)
mask_CN_4plot <- mask_CN[mask_CN$Chromosome=="chr3" & mask_CN$End>=9970000 & mask_CN$Start <=10130000,]
mask_CN_4plot$Start[1]<- 9970000
mask_CN_4plot$End[nrow(mask_CN_4plot)]<- 10130000
mask_FR_4plot <- mask_FR[mask_FR$Chromosome=="chr3" & mask_FR$End>=9970000 & mask_FR$Start <=10130000,]
mask_FR_4plot$Start[1]<- 9970000
mask_FR_4plot$End[nrow(mask_FR_4plot)]<- 10130000

p1<- ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_line(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3" & Set3.1_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std),size=2, alpha=1, colour = set_color[4])+
  geom_line(data = Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr3" & Set3.3_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std),size=2, alpha=1, colour = set_color[6])+
  # geom_line(data = Set3.2_xpclr_chrall[Set3.2_xpclr_chrall$Chr=="chr3" & Set3.2_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std),size=1.2, alpha=1)+
  # geom_line(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3" & SetAll_xpclr_chrall$PhysPos<=10130000,],aes(x=PhysPos, y=X3PCLR.B.std, color=set),size=1.2, alpha=1)+
  # geom_line()+
  scale_y_continuous(limits = c(-2.2,14), expand = c(0.1,0), breaks=seq(-2,14, by=2))+
  scale_x_continuous(limits=c(9970000, 10130000), expand = c(0.1,0), breaks = seq(9970000, 10130000, by=40000), labels=seq(9.97, 10.13, by=0.04))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  #unknown gene
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[4], y=c3.1B, yend=c3.1B, x=9970000, xend=10130000, linetype = 3, size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[6], y=c3.3B, yend=c3.3B, x=9970000, xend=10130000, linetype = 3, size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=9971744, xend=9975491,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=9993856, xend=9995638,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=11, yend=11, x=10008126, xend=10011788,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=1.7, yend=1.7, x=10012333, xend=10013691,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-0.2, yend=-0.2, x=10014619, xend=10016741,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=4, yend=4, x=10018411, xend=10017244,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-0.3, yend=-0.3, x=10022010, xend=10022120,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #unknown gene
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10047205, xend=10048334,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #trasnposon
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=0.5, yend=0.5, x=10075217, xend=10076201,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #NAC055
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10077093, xend=10078587,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+
  #NAC056
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",],y=-1, yend=-1, x=10106697, xend=10108574,color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=2)+

  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  # scale_color_manual(name="",breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],y=-2.2 - (14+2.2)*0.1, yend=-2.2 - (14+2.2)*0.1, x=9970000, xend=10130000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr8",],y=-2, yend=14, x=9970000+ (9970000 -10130000)*0.1, xend=9970000+ (9970000 -10130000)*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  geom_rect(data = mask_CN_4plot, aes(xmin=Start, xmax=End, ymin=-2.2, ymax=-1.6, fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  geom_rect(data = mask_FR_4plot, aes(xmin=Start, xmax=End, ymin=13.4, ymax=14,fill=type.0.single.copy.1.multi.copy),size=2, color="black", linewidth=0.05)+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=45), axis.title=element_text(size=45), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=45), legend.title = element_text(size=45), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(1,1,0,1), "cm"),  plot.title = element_text(size=45, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  annotate(geom = "text",x =  9971744+(9975491 -9971744)/2, y = -0.3, label=expression(italic("CKL8")), size=8, parse=T)+
  annotate(geom = "text",x =  9993856+(9995638 -9993856)/2, y = -0.3, label=expression(italic("AA3G18800")), size=8, parse=T)+
  annotate(geom = "text",x =  10008126+(10011788 -10008126)/2, y = 11.7, label=expression(italic("WDR44")), size=8, parse=T)+
  annotate(geom = "text",x =  10012333+(10013691 -10012333)/2, y =  2.4, label=expression(italic("AA3G18820")), size=8, parse=T)+
  annotate(geom = "text",x =  10014619+(10016741 -10014619)/2, y =  0.7, label=expression(italic("BRIX1-1")), size=8, parse=T)+
  annotate(geom = "text",x =  10017244+(10018411 -10017244)/2, y = 4.7, label=expression(italic("TSJT1")), size=8, parse=T)+
  annotate(geom = "text",x =  10022010+(10022120 -10022010)/2, y = -1, label=expression(italic("AA3G18850")), size=8, parse=T)+
  annotate(geom = "text",x =  10047205+(10048334 -10047205)/2, y = -0.3, label=expression(italic("AA3G18860")), size=8, parse=T)+
  annotate(geom = "text",x =  10075217+(10076201 -10075217)/2, y = 1.2, label=expression(italic("hAT-Ac")), size=8, parse=T)+
  annotate(geom = "text",x =  10077093+(10078587 -10077093)/2, y = -0.3, label=expression(italic("NAC055")), size=8, parse=T)+
  annotate(geom = "text",x =  10106697+(10108574 -10106697)/2, y = -0.3, label=expression(italic("NAC056")), size=8, parse=T)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))


p2<-ggplot()+
  # geom_point(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr4",],aes(x=PhysPos, y=X3PCLR.Anc.std, color=set))+
  geom_point(data = Set3.1_xpclr_chrall[Set3.1_xpclr_chrall$Chr=="chr3" ,],aes(x=PhysPos, y=X3PCLR.B.std),size=5, alpha=0.5, shape=21, color="black",fill = set_color[4])+
  geom_point(data = Set3.3_xpclr_chrall[Set3.3_xpclr_chrall$Chr=="chr3" ,],aes(x=PhysPos, y=X3PCLR.B.std),size=5, alpha=0.5, shape=21, color="black", fill = set_color[6])+
  # geom_line()+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[4], y=c3.1B, yend=c3.1B, x=0, xend=30421000, linetype = 3, size=2)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr1",], color=set_color[6], y=c3.3B, yend=c3.3B, x=0, xend=30421000, linetype = 3, size=2)+
  scale_y_continuous(limits = c(-2.2,14), expand = c(0.1,0), breaks=seq(-2,14, by=2))+
  scale_x_continuous(limits=c(0, 30421000), expand = c(0.1,0),breaks = seq(0, 30000000, by=5e6), labels=seq(0, 30, by=5))+
  # scale_color_gradient(low="blue3", high="brown2", breaks=c(0,1))+
  # scale_fill_manual(breaks = c("Set1", "Set2", "Set3", "Set4"), values = colors_sets)+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2.2 - (14+2.2)*0.1, yend=-2.2 - (14+2.2)*0.1, x=0, xend=30000000, color="black")+
  geom_segment(data = SetAll_xpclr_chrall[SetAll_xpclr_chrall$Chr=="chr3",],y=-2, yend=14, x=-30421000*0.1, xend=-30421000*0.1, color="black")+
  labs(x= "Position [Mbp]", y="CLR-score", title="")+
  scale_fill_manual(name="",breaks = c("0","1"), values=c("blue3", "brown2"), labels=c("single-copy", "multicopy"))+
  theme(panel.background = element_blank(), axis.text = element_text(size=45), axis.title=element_text(size=45), axis.ticks.length = unit(0.5, "cm"), legend.text = element_text(size=45), legend.title = element_text(size=45), legend.key.size = unit(1,"cm"), legend.position = "none", legend.key = element_rect(fill = "white",linetype = "blank"), legend.key.height = unit(1.5, "cm"), legend.key.width = unit(1.5, "cm"), plot.margin = unit(c(0,1,1,1), "cm"),  plot.title = element_text(size=45, hjust=0.5))+
  # annotate(geom = "text",x = 4.469e7, y =5.9, label="FR", size=12 , hjust=-0.5)+
  # annotate(geom = "text",x = 4.469e7, y =-2.1, label="CAN", size=12 , hjust=-0.4)+
  guides(color = guide_legend(override.aes = list(linewidth = 2)))



library(cowplot)
plot_grid(plotlist = list(p_SC5D, p1, p_chr6_S3.1_S2.1,p2 ), align = "hv")

