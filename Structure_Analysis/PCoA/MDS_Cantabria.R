library(ggplot2)
library(ggsci)
library(viridis)
library(pals)

safe_colorblind_gradient <- c("white", "#00009F", "#0000AF", "#0000BF", "#0000CF", "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", "#AF0000", "#9F0000", "#8F0000", "#800000")


safe_colorblind_palette <-c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                            "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "black", "orange", "purple", "pink","grey")

id2sample<-read.table("Final_Cantabrians.txt", stringsAsFactors = FALSE, header=F)


c15<-c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")
filled_shapes <- c(21,22,23,24,25,21,22,23,24,25,21,22,23,24)

# Load necessary library
library(MASS)  # For isoMDS to calculate stress

# Read the upper triangular distance matrix
x <- scan("CantabrianParaOutSNPable10PercentMissingness_pruned_clustering.mdist")
dims <- floor(sqrt(length(x) * 2))+1

m <- matrix(NA, dims, dims)
m[upper.tri(m, diag = F)] <- x
m <- t(m)

colnames(m) <- id2sample$V1
rownames(m) <- id2sample$V1



m<-as.dist(m)
# Ensure diagonal is zero (if it was not originally included)

# Perform Classical MDS
mds_result <- cmdscale(m, k = 2, eig = TRUE)
explained_variance <- mds_result$eig/sum(mds_result$eig[mds_result$eig>0])
sum(explained_variance[explained_variance>0])
# Extract coordinates
mds_coords <- mds_result$points
mds_coords<- as.data.frame(mds_coords)
mds_coords$population <- id2sample$V6[match(id2sample$V1, rownames(mds_coords))]
# Plot MDS result
ggplot(data=mds_coords, aes(x=V1,y=V2,fill=population, shape=population))+
  geom_point(size=7, color="black", stroke=2)+
  scale_fill_manual(values = c15)+
  scale_x_continuous(limits=c(-0.0475, 0.025), breaks=seq(-0.04, 0.02, by=0.02 ), expand=c(0,0))+
  scale_y_continuous(limits=c(-0.025, 0.03), breaks=seq(-0.02, 0.02, by=0.02 ), expand=c(0,0))+
  scale_shape_manual(values = filled_shapes)+
  geom_segment(y=-0.025, yend=-0.025, x=-0.04, xend=0.02)+
  geom_segment(y=-0.02, yend=0.02, x=-0.0475, xend=-0.0475)+
  labs(x="Coordinate 1 (14.1%)",y="Coordinate 2 (9.4%)", fill="Population", shape="Population")+
  theme(axis.title = element_text(size=40), legend.key.size = unit(1,"line"), axis.ticks.length = unit(0.5,"cm"),axis.text = element_text(size=40), legend.position = 0, panel.background = element_blank())
 
