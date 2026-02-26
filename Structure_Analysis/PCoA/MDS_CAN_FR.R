library(ggplot2)
library(ggsci)
library(viridis)
library(pals)

safe_colorblind_gradient <- c("white", "#00009F", "#0000AF", "#0000BF", "#0000CF", "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", "#AF0000", "#9F0000", "#8F0000", "#800000")

length(safe_colorblind_gradient)

safe_colorblind_palette <-c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                            "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "black", "orange", "purple", "pink","grey")

id2sample<-read.table("SampleList_FR_CAN_all_clean.txt", stringsAsFactors = FALSE, header=F)

mds<-read.table("/Structure_CAN_FR/plink.mds", stringsAsFactors = FALSE, header = TRUE)

mds$FID<-id2sample$V6[match(id2sample$V1,mds$IID)]
mds$FID <- factor(mds$FID,levels =  unique(mds$FID))
mds$REG <- "FR"
mds$REG[grep(x = mds$FID, pattern  = "^ES")] <- "ES"

c15<-c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")
filled_shapes <- c(21,22,23,24,25,21,22,23,24,25,21,22,23,24)
inner_shapes <- c(15,15,15,15,15,16,16,16,16,16,17,17,17,17,17)
inner_shapes <- c(21,21,21,21,21,22,22,22,22,22,23,23,23,23,23)
filled_shapes <- c(filled_shapes, filled_shapes, filled_shapes)
inner_shapes <-c(inner_shapes,inner_shapes,inner_shapes)
ggplot(data=mds, aes(x=C1,y=C2,fill=REG, shape=FID))+
  geom_point(size=7)+
  # scale_fill_manual(values = c15)+
  scale_shape_manual(values = filled_shapes)+
  labs(x="PCO 2",y="PCO 1", fill="Population", shape="Population")+
  theme(axis.title = element_text(size=15), legend.title = element_text(size =15), legend.key.size = unit(1,"line"),axis.text = element_text(size=15))+
  coord_fixed(ratio = 1)


  # Load necessary library
library(MASS)  # For isoMDS to calculate stress

# Read the upper triangular distance matrix
x <- scan("/home/btjeng/Data/Spanish_adaptation/Structure_CAN_FR/CAN_FR_pruned_clustering.mdist")
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
mds_coords
# Plot MDS result
mds_coords$REG<- "FR"
mds_coords$REG[grep(x = mds_coords$population, pattern  = "^ES")] <- "ES"

esp.col = rgb(255/255,164/255,5/255, 0.6)
scand.col= rgb(0,117/255,220/255, 0.6)
ceu.col=rgb(102/255,221/255,170/255, 0.6)
oc.col = rgb(255/255, 99/255, 71/255, 0.6)
length(unique(mds_coords$population))
length(grep(unique(mds_coords$population), pattern="FR"))
mds_coords<- rbind(mds_coords, mds_coords)
mds_coords$population <- factor(mds_coords$population, levels = unique(mds_coords$population))
p1<-ggplot(data=mds_coords, aes(x=V1,y=V2,fill=population, shape=population))+
  geom_point(size=7, color="black", stroke=2)+
  scale_fill_manual(values = c(rep(ceu.col,15 ),rep(esp.col,14)))+
  scale_x_continuous(limits=c(-0.04, 0.04), breaks=seq(-0.04, 0.04, by=0.02 ), expand=c(0.1,0))+
  scale_y_continuous(limits=c(-0.04, 0.02), breaks=seq(-0.04, 0.02, by=0.02 ), expand=c(0.1,0))+
  scale_shape_manual(values = filled_shapes)+
  geom_segment(y=-0.04, yend=0.02, x=-0.048, xend=-0.048)+
  geom_segment(y=-0.046, yend=-0.046, x=-0.04, xend=0.04)+
  labs(x="Coordinate 1 (42.3%)",y="Coordinate 2 (7.4%)", fill="Population", shape="Population")+
  theme(axis.title = element_text(size=40), legend.key.size = unit(1,"line"), axis.ticks.length = unit(0.5,"cm"),axis.text = element_text(size=40), legend.position = "right", panel.background = element_blank(),legend.text = element_text(size=30), legend.title = element_blank())


p2<-ggplot(data=mds_coords, aes(x=V1,y=V2,fill=population, shape=population))+
  geom_point(size=2, color="black", stroke=2)+
  scale_fill_manual(values = c(rep(ceu.col,15 ),rep(esp.col,14)))+
  scale_x_continuous(limits=c(-0.04, 0.04), breaks=seq(-0.04, 0.04, by=0.02 ), expand=c(0.1,0))+
  scale_y_continuous(limits=c(-0.04, 0.02), breaks=seq(-0.04, 0.02, by=0.02 ), expand=c(0.1,0))+
  scale_shape_manual(values = inner_shapes)+
  geom_segment(y=-0.04, yend=0.02, x=-0.048, xend=-0.048)+
  geom_segment(y=-0.046, yend=-0.046, x=-0.04, xend=0.04)+
  labs(x="Coordinate 1 (42.3%)",y="Coordinate 2 (7.4%)", fill="Population", shape="Population")+
  theme(axis.title = element_text(size=40), legend.key.size = unit(1,"line"), axis.ticks.length = unit(0.5,"cm"),axis.text = element_text(size=40), legend.position = "right", panel.background = element_blank(),legend.text = element_text(size=30), legend.title = element_blank())

library(ggplot2)
# ensure factor order matches your shape vectors
pops <- levels(mds_coords$population)

inner_map <- setNames(inner_shapes[seq_along(pops)], pops)

mds_coords$inner_pch <- inner_map[as.character(mds_coords$population)]
p <- p1 +
  geom_point(
    data = mds_coords,
    aes(x = V1, y = V2, fill = population),
    pch = mds_coords$inner_pch,  # <- literal pch vector
    size = 2,
    colour = "black",
    stroke = 1.5,
    inherit.aes = FALSE,
    show.legend = FALSE          # legend stays from p1
  )

p















ggplot(mds_coords, aes(V1, V2, fill = population, shape = population)) +
  geom_point(size = 7, color = "black", stroke = 2) +
  geom_point(size = 2, color = "black", stroke = 2) +
  scale_fill_manual(values = c(rep(ceu.col,15), rep(esp.col,14))) +
  scale_shape_manual(values = filled_shapes) +   # outer shapes (only one shape scale allowed)
  scale_x_continuous(limits=c(-0.04,0.04), breaks=seq(-0.04,0.04,0.02), expand=c(0.1,0)) +
  scale_y_continuous(limits=c(-0.04,0.02), breaks=seq(-0.04,0.02,0.02), expand=c(0.1,0)) +
  geom_segment(y=-0.04, yend=0.02, x=-0.048, xend=-0.048) +
  geom_segment(y=-0.046, yend=-0.046, x=-0.04, xend=0.04) +
  labs(x="Coordinate 1 (42.3%)", y="Coordinate 2 (7.4%)", fill="Population", shape="Population") +
  theme(axis.title = element_text(size=40),
        legend.key.size = unit(1,"line"),
        axis.ticks.length = unit(0.5,"cm"),
        axis.text = element_text(size=40),
        legend.position = "right",
        panel.background = element_blank(),
        legend.text = element_text(size=30),
        legend.title = element_blank())






#################


# 1) make sure population is a factor with stable order
mds_coords$population <- factor(mds_coords$population)

# 2) one letter per population level (NOT per point)
pop_levels  <- levels(mds_coords$population)
pop_letters <- c(letters[1:15],letters[1:14])   # or any vector you want
names(pop_letters) <- pop_levels
mds_coords$pop_letter <- pop_letters[as.character(mds_coords$population)]

ggplot(mds_coords, aes(V1, V2, fill = population, shape = population)) +
  geom_point(size = 7, color = "black", stroke = 2) +
  geom_text(aes(label = pop_letter), show.legend = FALSE) +  # <- per-point labels
  scale_fill_manual(values = c(rep(ceu.col,15), rep(esp.col,14))) +
  scale_shape_manual(values = filled_shapes) +
  guides(
    fill  = guide_legend(override.aes = list(label = pop_letters)),
    shape = guide_legend(override.aes = list(label = pop_letters))
  ) +
  scale_x_continuous(limits=c(-0.04,0.04), breaks=seq(-0.04,0.04,0.02), expand=c(0.1,0)) +
  scale_y_continuous(limits=c(-0.04,0.02), breaks=seq(-0.04,0.02,0.02), expand=c(0.1,0)) +
  geom_segment(y=-0.04, yend=0.02, x=-0.048, xend=-0.048) +
  geom_segment(y=-0.046, yend=-0.046, x=-0.04, xend=0.04) +
  labs(x="Coordinate 1 (42.3%)", y="Coordinate 2 (7.4%)", fill="Population", shape="Population") +
  theme(axis.title = element_text(size=40),
        legend.key.size = unit(1,"line"),
        axis.ticks.length = unit(0.5,"cm"),
        axis.text = element_text(size=40),
        legend.position = "right",
        panel.background = element_blank(),
        legend.text = element_text(size=30),
        legend.title = element_blank())




library(ggplot2)
library(grid)

# custom legend key: draw shape-21 point + centered label
draw_key_point_label <- function(data, params, size) {
  pt <- draw_key_point(data, params, size)
  tg <- textGrob(
    label = data$label %||% "",
    x = 0.5, y = 0.5,
    gp = gpar(col = data$colour %||% "black", fontsize = (data$size %||% 3.5) * 2)
  )
  grobTree(pt, tg)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# one letter per population level
mds_coords$population <- factor(mds_coords$population)
pop_levels  <- levels(mds_coords$population)
pop_letters <- setNames(pop_letters <- c(letters[1:15],letters[1:14]), pop_levels)
# pop_letters <- setNames(pop_letters <- c(rep("+", 5), rep("#", 5), rep("%", 5), rep("+", 5), rep("#", 5), rep("%", 4)), pop_levels)

# per-row label for plotting
mds_coords$pop_letter <- pop_letters[as.character(mds_coords$population)]

ggplot(mds_coords, aes(V1, V2, fill = population, shape = population)) +
  geom_point(
    # aes(label = pop_letter),                 # label aesthetic ONLY for the legend glyph
    size = 20, colour = "black", stroke = 1, alpha=0.5,
    key_glyph = draw_key_point_label
  ) +
  geom_text(aes(label = pop_letter), size=12, show.legend = FALSE, alpha=0.3) +
  scale_fill_manual(values = c(rep(ceu.col,15), rep(esp.col,14))) +
  scale_shape_manual(values = rep(21, 29)) +
  guides(
    shape = guide_legend(override.aes = list(size=20,label = unname(pop_letters))),
    fill  = guide_legend(override.aes = list(size=20,label = unname(pop_letters)))
  )+
  scale_x_continuous(limits=c(-0.036,0.05), breaks=seq(-0.03,0.03,0.01), expand=c(0,0)) +
  scale_y_continuous(limits=c(-0.04,0.01), breaks=seq(-0.04,0.01,0.01), expand=c(0.1,0)) +
  geom_segment(y=-0.04, yend=0.01, x=-0.036, xend=-0.036) +
  geom_segment(y=-0.045, yend=-0.045, x=-0.03, xend=0.03) +
  labs(x="Coordinate 1 (42.3%)", y="Coordinate 2 (7.4%)", fill="Population", shape="Population") +
  theme(axis.title = element_text(size=40),
        legend.key.size = unit(2,"line"),
        axis.ticks.length = unit(0.5,"cm"),
        axis.text = element_text(size=40),
        legend.position = c(0.9,0.5),
        panel.background = element_blank(),
        legend.text = element_text(size=30),
        legend.title = element_blank())


###

WCAN1<- c("ES01", "ES02", "ES03", "ES06","ES10")
WCAN2<- c("ES04", "ES05", "ES23")

ES17 <- "ES17"
FR01 <- "FR01"
CECAN <- c("ES24", "ES25", "ES08", "ES09")

mds_coords$SUBREG <- "FR"
mds_coords$SUBREG[mds_coords$population=="FR01"] <- "FR01"
mds_coords$SUBREG[mds_coords$population=="ES17"] <- "ES17"
mds_coords$SUBREG[mds_coords$population %in% CECAN] <-"CECAN"
mds_coords$SUBREG[mds_coords$population %in% c(WCAN1, WCAN2,"ES15")] <-"WCAN"

mds_coords$SUBREG <- factor(mds_coords$SUBREG, levels = c("FR01", "FR", "CECAN", "ES17", "WCAN"))
levels(mds_coords$SUBREG)

mds_coords_plot <- rbind(
  subset(mds_coords, SUBREG != "FR01"),
  subset(mds_coords, SUBREG == "FR01")
)

ggplot(mds_coords_plot, aes(V1, V2, fill = SUBREG, shape = SUBREG)) +
  geom_point(
    # aes(label = pop_letter),                 # label aesthetic ONLY for the legend glyph
    size = 15, colour = "black", stroke = 1, alpha=0.5,
  ) +
  scale_fill_manual(values = c(rep(ceu.col,2), rep(esp.col,3))) +
  scale_shape_manual(values = c(21,22,23,24,25)) +
  scale_x_continuous(limits=c(-0.036,0.05), breaks=seq(-0.03,0.03,0.01), expand=c(0,0)) +
  scale_y_continuous(limits=c(-0.04,0.01), breaks=seq(-0.04,0.01,0.01), expand=c(0.1,0)) +
  geom_segment(y=-0.04, yend=0.01, x=-0.036, xend=-0.036) +
  geom_segment(y=-0.045, yend=-0.045, x=-0.03, xend=0.03) +
  guides(
    shape = guide_legend(override.aes = list(alpha = 1), size=20),
    fill  = guide_legend(override.aes = list(alpha = 1),size=20)
  )+
  labs(x="Coordinate 1 (42.3%)", y="Coordinate 2 (7.4%)", fill="Population", shape="Population") +
  theme(axis.title = element_text(size=40),
        legend.key.size = unit(3,"line"),
        axis.ticks.length = unit(0.5,"cm"),
        axis.text = element_text(size=40),
        legend.position = c(0.9,0.5),
        panel.background = element_blank(),
        legend.text = element_text(size=40),
        legend.title = element_blank(), legend.key = element_rect(fill="white"))
