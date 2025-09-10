

flowering_data <- read.table("diedBefore_replicates_removed_excluding_FR04-264_FR04_258_paper_finalDataset.csv", header = T)

sample_map_Cantabria <-read.table("1000Genomes_final_Cantabrians.txt")
sample_map_FR <-read.table("SampleList_FR_all_clean.txt")

sample_map_Cantabria$genotype <- sapply(X = sample_map_Cantabria$V5, FUN  = function(x){paste(strsplit(x, split = "-")[[1]][1:2], collapse  = "-")})
sample_map_FR$genotype <- sapply(X = sample_map_FR$V5, FUN  = function(x){paste(strsplit(x, split = "-")[[1]][1:2], collapse  = "-")})
sample_map_all <- rbind(sample_map_Cantabria, sample_map_FR)

# sample_map_FR$genotype[sample_map_FR$genotype %in% flowering_data$genotype]
flowering_data$ID<-NA
for(i in 1:nrow(flowering_data)){
  flowering_data$ID[i] <- sample_map_all$V1[sample_map_all$genotype==flowering_data$genotype[i]]
}


write.table(flowering_data,file = "flowering_novern_seqid.txt", col.names = T, row.names = F, sep = "\t", quote = F)


mut <- read.table("FRL1_nonsyn_minGQ30DP5.tsv", header = F, sep = "\t")
colnames(mut) <- mut[1,]
mut <-mut[-1,]

mut$POS
flowering_data$mut_12618926 <- NA
flowering_data$mut_12618944 <- NA
flowering_data$mut_12618962 <- NA
flowering_data$mut_12618991 <- NA
flowering_data$mut_12618912 <- NA
i<-1
for(i in 1:nrow(flowering_data)){
  flowering_data[i,6:11] <- mut[,colnames(mut)==flowering_data$ID[i]]
}

sum(colnames(mut)==flowering_data$ID[i])

flowering_data


mds <- read.table(file = "plink_CanFR_tab.mds", sep="\t", header = T)
#
flowering_data$C1<- NA
flowering_data$C2<- NA
flowering_data$C3<- NA
flowering_data$C4<- NA

i<-1
for(i in 1:nrow(flowering_data)){
  flowering_data[i,12:15] <- mds[mds$FID==flowering_data$ID[i],c(4:7)]
}
flowering_data$DAS[flowering_data$DAS=="no_flower"] <-250
flowering_data$DAS <- as.numeric(flowering_data$DAS)

flowering_data[,6:11]

fit <- lm(data = flowering_data, formula = DAS~ mut_12618926)
plot(fit, which = 1)   # residuals vs fitted
plot(fit, which = 3)   # residuals vs fitted
summary(fit)
qqnorm(resid(fit)); qqline(resid(fit))

  library(lme4)
# install.packages("lmerTest")
library(lmerTest)

fit <- lmer(DAS ~ mut_12618926 + (1|genotype), data=flowering_data)
summary(fit)

fit1 <- lmer(DAS ~ mut_12618926 + (1|genotype), data = flowering_data)
fit0 <- lmer(DAS ~ 1 + (1|genotype), data = flowering_data)

anova(fit0, fit1)  # chi-sq test for fixed effect

library(car)
Anova(fit1, type=3)

plot(fitted(fit1), resid(fit1))
abline(h=0, col="red")
qqnorm(resid(fit1)); qqline(resid(fit1))

fit1 <- lmer(DAS ~ 1 + mut_12618926  +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)
fit0 <- lmer(DAS ~ 1   +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)


anova(fit0, fit1)  # chi-sq test for fixed effect

library(car)
Anova(fit1, type=3)

plot(fitted(fit1), resid(fit1))
abline(h=0, col="red")
qqnorm(resid(fit1)); qqline(resid(fit1))




##
fit1 <- lmer(DAS ~ 1 + mut_12618944  +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)
fit0 <- lmer(DAS ~ 1   +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)


anova(fit0, fit1)  # chi-sq test for fixed effect

fit1 <- lmer(DAS ~ 1 + mut_12619056  +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)
fit0 <- lmer(DAS ~ 1   +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)


anova(fit0, fit1)  # chi-sq test for fixed effect


fit1 <- lmer(DAS ~ 1 + mut_12618962  +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)
fit0 <- lmer(DAS ~ 1   +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)


anova(fit0, fit1)  # chi-sq test for fixed effect


fit1 <- lmer(DAS ~ 1 + mut_12618991  +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)
fit0 <- lmer(DAS ~ 1   +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)


anova(fit0, fit1)  # chi-sq test for fixed effect


fit1 <- lmer(DAS ~ 1 + mut_12618912  +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)
fit0 <- lmer(DAS ~ 1   +C1 + C2 + C3 + C4 + (1|genotype), data = flowering_data)


######
# fit1 <- lmer(DAS ~ 1   +(1|C1) +(1|C2)+(1|C3)+(1|C4)+ (1|genotype), data = flowering_data)
# fit0 <- lmer(DAS ~ 1   +(1|C1) +(1|C2)+(1|C3)+(1|C4)+ (1|genotype), data = flowering_data)
#

anova(fit0, fit1)  # chi-sq test for fixed effect

