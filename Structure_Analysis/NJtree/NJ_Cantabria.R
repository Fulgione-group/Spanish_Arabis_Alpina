library(ape)

x <- scan("CantabrianParaOutSNPable10PercentMissingness_pruned_clustering.mdist")
dims <- floor(sqrt(length(x) * 2))+1

m <- matrix(NA, dims, dims)
m[upper.tri(m, diag = F)] <- x
m <- t(m)


# Get the ordered sample names
# Take out 1 column of the <prefix>.mdist.id file to make samples_<prefix>.txt file

samples <- read.table("CantabrianParaOutSNPable10PercentMissingness_pruned_clustering.mdist.id")
samples<- samples$V1

sampleAnn<- read.table("1000Genomes_final_Iberians.txt", header=F, sep = "\t")
head(sampleAnn)
pops<-sampleAnn$V6[match(
  samples,
  sampleAnn$V1
)]

# Add sample names to the matrix
# This is probably not strictly necessary but I needed it for the metadata

dim(m)
m[1,]
colnames(m) <- samples
rownames(m) <- samples

dist <- as.dist(m)

# Make the tree
xx <- as.dist(as(dist, "matrix"))
trw <- nj(xx)

## If you do not need colours and only sample names, you can skip to "End skip" ##
# Prepare plotting
id <- as.character(trw$tip.label)
ord <- seq(1, length(id))
# I convert the list of tip labels (which should just be the sample names) into a data frame

tlab <- as.data.frame(cbind(id, ord))

# Because tlab is now a data frame, I can merge it with my metadata file here
library(RColorBrewer2)
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

c15<-c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")
filled_shapes <- c(21,22,23,24,25,21,22,23,24,25,21,22,23,24)

cdf<- data.frame(pop=unique(pops),colour= c15[1:14], shape =filled_shapes)

labs <- cbind(tlab,pops)
labs$colour <- cdf$colour[match(labs$pops,cdf$pop)]
labs$shape<- cdf$shape[match(labs$pops,cdf$pop)]
iro <- as.character(labs$colour)


# You could alternatively not use a metadata file and make a vector of hex colour code (as character) of length n(samples)

## End skip ##


# Plot
# Do not specify tip.color if you don't need them
plot(trw, "unrooted", main="", cex=1, tip.color = iro, no.margin = T)
tiplabels(pch=labs$shape, col=iro, cex=3)

legend(x = 0.14, y = 0.1, legend=cdf$pop,
       col=cdf$colour, cex=0.5, lty = 1)

# If you do not want coloured sample names like me, you could do this to have colours only:
plot(trw, "unrooted", main="", cex=0.25, tip.color = iro, no.margin = T, show.tip.label = F)
tiplabels(pch=labs$shape, col=iro, cex=1)




plot(trw, "unrooted", main="", cex=0.25, tip.color = iro, no.margin = T,show.tip.label=F)
tiplabels(pch=labs$shape, col="black", cex=2, bg=iro, lwd=2)
legend(x = 0.09, y = 0.075, legend=cdf$pop,
       ,pt.bg=cdf$colour, col="black",pch=cdf$shape, cex=0.75, pt.cex = 2, text.font = 2,  pt.lwd=2)
