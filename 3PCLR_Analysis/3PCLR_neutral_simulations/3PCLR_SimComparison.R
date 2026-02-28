#set3.1



#start with set4

#Load simulated 3P-CLR
sim_wcan1_wcan2_CECAN<-read.table("simulated_repall.vcf.WCAN1.WCAN2.CECAN.chr1.1PerCentMaf.3PCLRoutput",sep = "\t", header = T)
#Load emirical 3P-CLR
Set3.1_xpclr_chrall<- read.table(Set3.1_xpclr_chrall_new.txt", sep = "\t", header=T)


# standardize simulated scores
sim_wcan1_wcan2_CECAN$X3PCLR.Anc.std <- (sim_wcan1_wcan2_CECAN$X3PCLR.Anc-mean(sim_wcan1_wcan2_CECAN$X3PCLR.Anc))/sd(sim_wcan1_wcan2_CECAN$X3PCLR.Anc)
sim_wcan1_wcan2_CECAN$X3PCLR.A.std <- (sim_wcan1_wcan2_CECAN$X3PCLR.A-mean(sim_wcan1_wcan2_CECAN$X3PCLR.A))/sd(sim_wcan1_wcan2_CECAN$X3PCLR.A)
sim_wcan1_wcan2_CECAN$X3PCLR.B.std <- (sim_wcan1_wcan2_CECAN$X3PCLR.B-mean(sim_wcan1_wcan2_CECAN$X3PCLR.B))/sd(sim_wcan1_wcan2_CECAN$X3PCLR.B)

#Compare quantiles (5%,1%,0.5%)

#WCAN1
hist(sim_wcan1_wcan2_CECAN$X3PCLR.A)


sort(Set3.1_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.05*length(Set3.1_xpclr_chrall$X3PCLR.B)]
sort(Set3.1_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.01*length(Set3.1_xpclr_chrall$X3PCLR.B)]
sort(Set3.1_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.B)]

sort(Set3.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.05*length(Set3.1_xpclr_chrall$X3PCLR.B.std)]
sort(Set3.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set3.1_xpclr_chrall$X3PCLR.B.std)]
sort(Set3.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.B.std)]

sort(sim_wcan1_wcan2_CECAN$X3PCLR.A, decreasing = T)[0.05*length(sim_wcan1_wcan2_CECAN$X3PCLR.A)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.A, decreasing = T)[0.01*length(sim_wcan1_wcan2_CECAN$X3PCLR.A)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.A, decreasing = T)[0.005*length(sim_wcan1_wcan2_CECAN$X3PCLR.A)]

sort(sim_wcan1_wcan2_CECAN$X3PCLR.A.std, decreasing = T)[0.05*length(sim_wcan1_wcan2_CECAN$X3PCLR.A.std)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.A.std, decreasing = T)[0.01*length(sim_wcan1_wcan2_CECAN$X3PCLR.A.std)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.A.std, decreasing = T)[0.005*length(sim_wcan1_wcan2_CECAN$X3PCLR.A.std)]

#WCAN2
hist(sim_wcan1_wcan2_CECAN$X3PCLR.B)

sort(Set3.1_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.05*length(Set3.1_xpclr_chrall$X3PCLR.A)]
sort(Set3.1_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.01*length(Set3.1_xpclr_chrall$X3PCLR.A)]
sort(Set3.1_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.A)]

sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.05*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]
sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]
sort(Set3.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.A.std)]

sort(sim_wcan1_wcan2_CECAN$X3PCLR.B, decreasing = T)[0.05*length(sim_wcan1_wcan2_CECAN$X3PCLR.B)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.B, decreasing = T)[0.01*length(sim_wcan1_wcan2_CECAN$X3PCLR.B)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.B, decreasing = T)[0.005*length(sim_wcan1_wcan2_CECAN$X3PCLR.B)]

sort(sim_wcan1_wcan2_CECAN$X3PCLR.B.std, decreasing = T)[0.05*length(sim_wcan1_wcan2_CECAN$X3PCLR.B.std)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.B.std, decreasing = T)[0.01*length(sim_wcan1_wcan2_CECAN$X3PCLR.B.std)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.B.std, decreasing = T)[0.005*length(sim_wcan1_wcan2_CECAN$X3PCLR.B.std)]

#ANC
hist(sim_wcan1_wcan2_CECAN$X3PCLR.Anc)

sort(Set3.1_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.05*length(Set3.1_xpclr_chrall$X3PCLR.Anc)]
sort(Set3.1_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.01*length(Set3.1_xpclr_chrall$X3PCLR.Anc)]
sort(Set3.1_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.Anc)]

sort(Set3.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.05*length(Set3.1_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set3.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set3.1_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set3.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set3.1_xpclr_chrall$X3PCLR.Anc.std)]

sort(sim_wcan1_wcan2_CECAN$X3PCLR.Anc, decreasing = T)[0.05*length(sim_wcan1_wcan2_CECAN$X3PCLR.Anc)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.Anc, decreasing = T)[0.01*length(sim_wcan1_wcan2_CECAN$X3PCLR.Anc)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.Anc, decreasing = T)[0.005*length(sim_wcan1_wcan2_CECAN$X3PCLR.Anc)]

sort(sim_wcan1_wcan2_CECAN$X3PCLR.Anc.std, decreasing = T)[0.05*length(sim_wcan1_wcan2_CECAN$X3PCLR.Anc.std)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.Anc.std, decreasing = T)[0.01*length(sim_wcan1_wcan2_CECAN$X3PCLR.Anc.std)]
sort(sim_wcan1_wcan2_CECAN$X3PCLR.Anc.std, decreasing = T)[0.005*length(sim_wcan1_wcan2_CECAN$X3PCLR.Anc.std)]




#Set4

sim_WCAN2_ES17_CECAN1<-read.table("simulated_repall.vcf.WCAN2.ES17.CECAN1.chr1.1PerCentMaf.3PCLRoutput", sep = "\t", header = T)
Set3.2_xpclr_chrall<- read.table(Set3.2_xpclr_chrall_new.txt", sep = "\t", header=T)

sim_WCAN2_ES17_CECAN1$X3PCLR.Anc.std <- (sim_WCAN2_ES17_CECAN1$X3PCLR.Anc-mean(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc))/sd(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc)
sim_WCAN2_ES17_CECAN1$X3PCLR.A.std <- (sim_WCAN2_ES17_CECAN1$X3PCLR.A-mean(sim_WCAN2_ES17_CECAN1$X3PCLR.A))/sd(sim_WCAN2_ES17_CECAN1$X3PCLR.A)
sim_WCAN2_ES17_CECAN1$X3PCLR.B.std <- (sim_WCAN2_ES17_CECAN1$X3PCLR.B-mean(sim_WCAN2_ES17_CECAN1$X3PCLR.B))/sd(sim_WCAN2_ES17_CECAN1$X3PCLR.B)

#WCAN2
hist(sim_WCAN2_ES17_CECAN1$X3PCLR.A)

sort(Set3.2_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.05*length(Set3.2_xpclr_chrall$X3PCLR.A)]
sort(Set3.2_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.01*length(Set3.2_xpclr_chrall$X3PCLR.A)]
sort(Set3.2_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.A)]

sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.05*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]
sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]
sort(Set3.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.A.std)]

sort(sim_WCAN2_ES17_CECAN1$X3PCLR.A, decreasing = T)[0.05*length(sim_WCAN2_ES17_CECAN1$X3PCLR.A)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.A, decreasing = T)[0.01*length(sim_WCAN2_ES17_CECAN1$X3PCLR.A)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.A, decreasing = T)[0.005*length(sim_WCAN2_ES17_CECAN1$X3PCLR.A)]

sort(sim_WCAN2_ES17_CECAN1$X3PCLR.A.std, decreasing = T)[0.05*length(sim_WCAN2_ES17_CECAN1$X3PCLR.A.std)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.A.std, decreasing = T)[0.01*length(sim_WCAN2_ES17_CECAN1$X3PCLR.A.std)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.A.std, decreasing = T)[0.005*length(sim_WCAN2_ES17_CECAN1$X3PCLR.A.std)]






#ES17
hist(sim_WCAN2_ES17_CECAN1$X3PCLR.B)

sort(Set3.2_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.05*length(Set3.2_xpclr_chrall$X3PCLR.B)]
sort(Set3.2_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.01*length(Set3.2_xpclr_chrall$X3PCLR.B)]
sort(Set3.2_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.B)]

sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.05*length(Set3.2_xpclr_chrall$X3PCLR.B.std)]
sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set3.2_xpclr_chrall$X3PCLR.B.std)]
sort(Set3.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.B.std)]

sort(sim_WCAN2_ES17_CECAN1$X3PCLR.B, decreasing = T)[0.05*length(sim_WCAN2_ES17_CECAN1$X3PCLR.B)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.B, decreasing = T)[0.01*length(sim_WCAN2_ES17_CECAN1$X3PCLR.B)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.B, decreasing = T)[0.005*length(sim_WCAN2_ES17_CECAN1$X3PCLR.B)]

sort(sim_WCAN2_ES17_CECAN1$X3PCLR.B.std, decreasing = T)[0.05*length(sim_WCAN2_ES17_CECAN1$X3PCLR.B.std)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.B.std, decreasing = T)[0.01*length(sim_WCAN2_ES17_CECAN1$X3PCLR.B.std)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.B.std, decreasing = T)[0.005*length(sim_WCAN2_ES17_CECAN1$X3PCLR.B.std)]






#ANC
hist(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc)

sort(Set3.2_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.05*length(Set3.2_xpclr_chrall$X3PCLR.Anc)]
sort(Set3.2_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.01*length(Set3.2_xpclr_chrall$X3PCLR.Anc)]
sort(Set3.2_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.Anc)]

sort(Set3.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.05*length(Set3.2_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set3.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set3.2_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set3.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set3.2_xpclr_chrall$X3PCLR.Anc.std)]

sort(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc, decreasing = T)[0.05*length(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc, decreasing = T)[0.01*length(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc, decreasing = T)[0.005*length(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc)]

sort(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc.std, decreasing = T)[0.05*length(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc.std)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc.std, decreasing = T)[0.01*length(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc.std)]
sort(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc.std, decreasing = T)[0.005*length(sim_WCAN2_ES17_CECAN1$X3PCLR.Anc.std)]



#Set6
sim_WCAN1_ES17_CECAN1<-read.table("~/Data/Spanish_adaptation//simulations_output/simulated_repall.vcf.WCAN1.ES17.CECAN1.chr1.1PerCentMaf.3PCLRoutput", sep = "\t", header = T)
Set3.3_xpclr_chrall<- read.table(Set3.3_xpclr_chrall_new.txt", sep = "\t", header=T)

sim_WCAN1_ES17_CECAN1$X3PCLR.Anc.std <- (sim_WCAN1_ES17_CECAN1$X3PCLR.Anc-mean(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc))/sd(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc)
sim_WCAN1_ES17_CECAN1$X3PCLR.A.std <- (sim_WCAN1_ES17_CECAN1$X3PCLR.A-mean(sim_WCAN1_ES17_CECAN1$X3PCLR.A))/sd(sim_WCAN1_ES17_CECAN1$X3PCLR.A)
sim_WCAN1_ES17_CECAN1$X3PCLR.B.std <- (sim_WCAN1_ES17_CECAN1$X3PCLR.B-mean(sim_WCAN1_ES17_CECAN1$X3PCLR.B))/sd(sim_WCAN1_ES17_CECAN1$X3PCLR.B)

#WCAN1
hist(sim_WCAN1_ES17_CECAN1$X3PCLR.A)

sort(Set3.3_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.05*length(Set3.3_xpclr_chrall$X3PCLR.A)]
sort(Set3.3_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.01*length(Set3.3_xpclr_chrall$X3PCLR.A)]
sort(Set3.3_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.005*length(Set3.3_xpclr_chrall$X3PCLR.A)]

sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.05*length(Set3.3_xpclr_chrall$X3PCLR.A.std)]
sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set3.3_xpclr_chrall$X3PCLR.A.std)]
sort(Set3.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set3.3_xpclr_chrall$X3PCLR.A.std)]

sort(sim_WCAN1_ES17_CECAN1$X3PCLR.A, decreasing = T)[0.05*length(sim_WCAN1_ES17_CECAN1$X3PCLR.A)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.A, decreasing = T)[0.01*length(sim_WCAN1_ES17_CECAN1$X3PCLR.A)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.A, decreasing = T)[0.005*length(sim_WCAN1_ES17_CECAN1$X3PCLR.A)]

sort(sim_WCAN1_ES17_CECAN1$X3PCLR.A.std, decreasing = T)[0.05*length(sim_WCAN1_ES17_CECAN1$X3PCLR.A.std)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.A.std, decreasing = T)[0.01*length(sim_WCAN1_ES17_CECAN1$X3PCLR.A.std)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.A.std, decreasing = T)[0.005*length(sim_WCAN1_ES17_CECAN1$X3PCLR.A.std)]






#ES17
hist(sim_WCAN1_ES17_CECAN1$X3PCLR.B)

sort(Set3.3_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.05*length(Set3.3_xpclr_chrall$X3PCLR.B)]
sort(Set3.3_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.01*length(Set3.3_xpclr_chrall$X3PCLR.B)]
sort(Set3.3_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.005*length(Set3.3_xpclr_chrall$X3PCLR.B)]

sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.05*length(Set3.3_xpclr_chrall$X3PCLR.B.std)]
sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set3.3_xpclr_chrall$X3PCLR.B.std)]
sort(Set3.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set3.3_xpclr_chrall$X3PCLR.B.std)]

sort(sim_WCAN1_ES17_CECAN1$X3PCLR.B, decreasing = T)[0.05*length(sim_WCAN1_ES17_CECAN1$X3PCLR.B)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.B, decreasing = T)[0.01*length(sim_WCAN1_ES17_CECAN1$X3PCLR.B)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.B, decreasing = T)[0.005*length(sim_WCAN1_ES17_CECAN1$X3PCLR.B)]

sort(sim_WCAN1_ES17_CECAN1$X3PCLR.B.std, decreasing = T)[0.05*length(sim_WCAN1_ES17_CECAN1$X3PCLR.B.std)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.B.std, decreasing = T)[0.01*length(sim_WCAN1_ES17_CECAN1$X3PCLR.B.std)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.B.std, decreasing = T)[0.005*length(sim_WCAN1_ES17_CECAN1$X3PCLR.B.std)]






#ANC
hist(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc)

sort(Set3.3_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.05*length(Set3.3_xpclr_chrall$X3PCLR.Anc)]
sort(Set3.3_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.01*length(Set3.3_xpclr_chrall$X3PCLR.Anc)]
sort(Set3.3_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.005*length(Set3.3_xpclr_chrall$X3PCLR.Anc)]

sort(Set3.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.05*length(Set3.3_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set3.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set3.3_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set3.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set3.3_xpclr_chrall$X3PCLR.Anc.std)]

sort(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc, decreasing = T)[0.05*length(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc, decreasing = T)[0.01*length(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc, decreasing = T)[0.005*length(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc)]

sort(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc.std, decreasing = T)[0.05*length(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc.std)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc.std, decreasing = T)[0.01*length(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc.std)]
sort(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc.std, decreasing = T)[0.005*length(sim_WCAN1_ES17_CECAN1$X3PCLR.Anc.std)]




#Set1

sim_wcan1_wcan2_fr<-read.table("~/Data/Spanish_adaptation//simulations_output/simulated_repall.vcf.WCAN1.WCAN2.FR.chr1.1PerCentMaf.3PCLRoutput", sep = "\t", header = T)
Set2.1_xpclr_chrall<- read.table(Set2.1_xpclr_chrall_new.txt", sep = "\t", header=T)


sim_wcan1_wcan2_fr$X3PCLR.Anc.std <- (sim_wcan1_wcan2_fr$X3PCLR.Anc-mean(sim_wcan1_wcan2_fr$X3PCLR.Anc))/sd(sim_wcan1_wcan2_fr$X3PCLR.Anc)
sim_wcan1_wcan2_fr$X3PCLR.A.std <- (sim_wcan1_wcan2_fr$X3PCLR.A-mean(sim_wcan1_wcan2_fr$X3PCLR.A))/sd(sim_wcan1_wcan2_fr$X3PCLR.A)
sim_wcan1_wcan2_fr$X3PCLR.B.std <- (sim_wcan1_wcan2_fr$X3PCLR.B-mean(sim_wcan1_wcan2_fr$X3PCLR.B))/sd(sim_wcan1_wcan2_fr$X3PCLR.B)


#WCAN1
hist(sim_wcan1_wcan2_fr$X3PCLR.A)

sort(Set2.1_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.05*length(Set2.1_xpclr_chrall$X3PCLR.B)]
sort(Set2.1_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.01*length(Set2.1_xpclr_chrall$X3PCLR.B)]
sort(Set2.1_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.B)]

sort(Set2.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.05*length(Set2.1_xpclr_chrall$X3PCLR.B.std)]
sort(Set2.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set2.1_xpclr_chrall$X3PCLR.B.std)]
sort(Set2.1_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.B.std)]

sort(sim_wcan1_wcan2_fr$X3PCLR.A, decreasing = T)[0.05*length(sim_wcan1_wcan2_CECAN$X3PCLR.A)]
sort(sim_wcan1_wcan2_fr$X3PCLR.A, decreasing = T)[0.01*length(sim_wcan1_wcan2_fr$X3PCLR.A)]
sort(sim_wcan1_wcan2_fr$X3PCLR.A, decreasing = T)[0.005*length(sim_wcan1_wcan2_fr$X3PCLR.A)]

sort(sim_wcan1_wcan2_fr$X3PCLR.A.std, decreasing = T)[0.05*length(sim_wcan1_wcan2_fr$X3PCLR.A.std)]
sort(sim_wcan1_wcan2_fr$X3PCLR.A.std, decreasing = T)[0.01*length(sim_wcan1_wcan2_fr$X3PCLR.A.std)]
sort(sim_wcan1_wcan2_fr$X3PCLR.A.std, decreasing = T)[0.005*length(sim_wcan1_wcan2_fr$X3PCLR.A.std)]

#WCAN2
hist(sim_wcan1_wcan2_fr$X3PCLR.B)

sort(Set2.1_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.05*length(Set2.1_xpclr_chrall$X3PCLR.A)]
sort(Set2.1_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.01*length(Set2.1_xpclr_chrall$X3PCLR.A)]
sort(Set2.1_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.A)]

sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.05*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]
sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]
sort(Set2.1_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.A.std)]

sort(sim_wcan1_wcan2_fr$X3PCLR.B, decreasing = T)[0.05*length(sim_wcan1_wcan2_fr$X3PCLR.B)]
sort(sim_wcan1_wcan2_fr$X3PCLR.B, decreasing = T)[0.01*length(sim_wcan1_wcan2_fr$X3PCLR.B)]
sort(sim_wcan1_wcan2_fr$X3PCLR.B, decreasing = T)[0.005*length(sim_wcan1_wcan2_fr$X3PCLR.B)]

sort(sim_wcan1_wcan2_fr$X3PCLR.B.std, decreasing = T)[0.05*length(sim_wcan1_wcan2_fr$X3PCLR.B.std)]
sort(sim_wcan1_wcan2_fr$X3PCLR.B.std, decreasing = T)[0.01*length(sim_wcan1_wcan2_fr$X3PCLR.B.std)]
sort(sim_wcan1_wcan2_fr$X3PCLR.B.std, decreasing = T)[0.005*length(sim_wcan1_wcan2_fr$X3PCLR.B.std)]

#ANC
hist(sim_wcan1_wcan2_fr$X3PCLR.Anc)

sort(Set2.1_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.05*length(Set2.1_xpclr_chrall$X3PCLR.Anc)]
sort(Set2.1_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.01*length(Set2.1_xpclr_chrall$X3PCLR.Anc)]
sort(Set2.1_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.Anc)]

sort(Set2.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.05*length(Set2.1_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set2.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set2.1_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set2.1_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set2.1_xpclr_chrall$X3PCLR.Anc.std)]

sort(sim_wcan1_wcan2_fr$X3PCLR.Anc, decreasing = T)[0.05*length(sim_wcan1_wcan2_fr$X3PCLR.Anc)]
sort(sim_wcan1_wcan2_fr$X3PCLR.Anc, decreasing = T)[0.01*length(sim_wcan1_wcan2_fr$X3PCLR.Anc)]
sort(sim_wcan1_wcan2_fr$X3PCLR.Anc, decreasing = T)[0.005*length(sim_wcan1_wcan2_fr$X3PCLR.Anc)]

sort(sim_wcan1_wcan2_fr$X3PCLR.Anc.std, decreasing = T)[0.05*length(sim_wcan1_wcan2_fr$X3PCLR.Anc.std)]
sort(sim_wcan1_wcan2_fr$X3PCLR.Anc.std, decreasing = T)[0.01*length(sim_wcan1_wcan2_fr$X3PCLR.Anc.std)]
sort(sim_wcan1_wcan2_fr$X3PCLR.Anc.std, decreasing = T)[0.005*length(sim_wcan1_wcan2_fr$X3PCLR.Anc.std)]





#Set2
sim_WCAN2_ES17_FR<-read.table("~/Data/Spanish_adaptation//simulations_output/simulated_repall.vcf.WCAN2.ES17.FR.chr1.1PerCentMaf.3PCLRoutput", sep = "\t", header = T)
Set2.2_xpclr_chrall<- read.table(Set2.2_xpclr_chrall_new.txt", sep = "\t", header=T)

sim_WCAN2_ES17_FR$X3PCLR.Anc.std <- (sim_WCAN2_ES17_FR$X3PCLR.Anc-mean(sim_WCAN2_ES17_FR$X3PCLR.Anc))/sd(sim_WCAN2_ES17_FR$X3PCLR.Anc)
sim_WCAN2_ES17_FR$X3PCLR.A.std <- (sim_WCAN2_ES17_FR$X3PCLR.A-mean(sim_WCAN2_ES17_FR$X3PCLR.A))/sd(sim_WCAN2_ES17_FR$X3PCLR.A)
sim_WCAN2_ES17_FR$X3PCLR.B.std <- (sim_WCAN2_ES17_FR$X3PCLR.B-mean(sim_WCAN2_ES17_FR$X3PCLR.B))/sd(sim_WCAN2_ES17_FR$X3PCLR.B)


#WCAN2
hist(sim_WCAN2_ES17_FR$X3PCLR.A)

sort(Set2.2_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.05*length(Set2.2_xpclr_chrall$X3PCLR.A)]
sort(Set2.2_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.01*length(Set2.2_xpclr_chrall$X3PCLR.A)]
sort(Set2.2_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.005*length(Set2.2_xpclr_chrall$X3PCLR.A)]

sort(Set2.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.05*length(Set2.2_xpclr_chrall$X3PCLR.A.std)]
sort(Set2.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set2.2_xpclr_chrall$X3PCLR.A.std)]
sort(Set2.2_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set2.2_xpclr_chrall$X3PCLR.A.std)]

sort(sim_WCAN2_ES17_FR$X3PCLR.A, decreasing = T)[0.05*length(sim_wcan1_wcan2_CECAN$X3PCLR.A)]
sort(sim_WCAN2_ES17_FR$X3PCLR.A, decreasing = T)[0.01*length(sim_WCAN2_ES17_FR$X3PCLR.A)]
sort(sim_WCAN2_ES17_FR$X3PCLR.A, decreasing = T)[0.005*length(sim_WCAN2_ES17_FR$X3PCLR.A)]

sort(sim_WCAN2_ES17_FR$X3PCLR.A.std, decreasing = T)[0.05*length(sim_WCAN2_ES17_FR$X3PCLR.A.std)]
sort(sim_WCAN2_ES17_FR$X3PCLR.A.std, decreasing = T)[0.01*length(sim_WCAN2_ES17_FR$X3PCLR.A.std)]
sort(sim_WCAN2_ES17_FR$X3PCLR.A.std, decreasing = T)[0.005*length(sim_WCAN2_ES17_FR$X3PCLR.A.std)]

#ES17
hist(sim_WCAN2_ES17_FR$X3PCLR.B)

sort(Set2.2_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.05*length(Set2.2_xpclr_chrall$X3PCLR.B)]
sort(Set2.2_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.01*length(Set2.2_xpclr_chrall$X3PCLR.B)]
sort(Set2.2_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.005*length(Set2.2_xpclr_chrall$X3PCLR.B)]

sort(Set2.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.05*length(Set2.2_xpclr_chrall$X3PCLR.B.std)]
sort(Set2.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set2.2_xpclr_chrall$X3PCLR.B.std)]
sort(Set2.2_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set2.2_xpclr_chrall$X3PCLR.B.std)]

sort(sim_WCAN2_ES17_FR$X3PCLR.B, decreasing = T)[0.05*length(sim_WCAN2_ES17_FR$X3PCLR.B)]
sort(sim_WCAN2_ES17_FR$X3PCLR.B, decreasing = T)[0.01*length(sim_WCAN2_ES17_FR$X3PCLR.B)]
sort(sim_WCAN2_ES17_FR$X3PCLR.B, decreasing = T)[0.005*length(sim_WCAN2_ES17_FR$X3PCLR.B)]

sort(sim_WCAN2_ES17_FR$X3PCLR.B.std, decreasing = T)[0.05*length(sim_WCAN2_ES17_FR$X3PCLR.B.std)]
sort(sim_WCAN2_ES17_FR$X3PCLR.B.std, decreasing = T)[0.01*length(sim_WCAN2_ES17_FR$X3PCLR.B.std)]
sort(sim_WCAN2_ES17_FR$X3PCLR.B.std, decreasing = T)[0.005*length(sim_WCAN2_ES17_FR$X3PCLR.B.std)]

#ANC
hist(sim_WCAN2_ES17_FR$X3PCLR.Anc)

sort(Set2.2_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.05*length(Set2.2_xpclr_chrall$X3PCLR.Anc)]
sort(Set2.2_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.01*length(Set2.2_xpclr_chrall$X3PCLR.Anc)]
sort(Set2.2_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.005*length(Set2.2_xpclr_chrall$X3PCLR.Anc)]

sort(Set2.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.05*length(Set2.2_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set2.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set2.2_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set2.2_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set2.2_xpclr_chrall$X3PCLR.Anc.std)]

sort(sim_WCAN2_ES17_FR$X3PCLR.Anc, decreasing = T)[0.05*length(sim_WCAN2_ES17_FR$X3PCLR.Anc)]
sort(sim_WCAN2_ES17_FR$X3PCLR.Anc, decreasing = T)[0.01*length(sim_WCAN2_ES17_FR$X3PCLR.Anc)]
sort(sim_WCAN2_ES17_FR$X3PCLR.Anc, decreasing = T)[0.005*length(sim_WCAN2_ES17_FR$X3PCLR.Anc)]

sort(sim_WCAN2_ES17_FR$X3PCLR.Anc.std, decreasing = T)[0.05*length(sim_WCAN2_ES17_FR$X3PCLR.Anc.std)]
sort(sim_WCAN2_ES17_FR$X3PCLR.Anc.std, decreasing = T)[0.01*length(sim_WCAN2_ES17_FR$X3PCLR.Anc.std)]
sort(sim_WCAN2_ES17_FR$X3PCLR.Anc.std, decreasing = T)[0.005*length(sim_WCAN2_ES17_FR$X3PCLR.Anc.std)]




#Set3
sim_WCAN1_ES17_FR<-read.table("~/Data/Spanish_adaptation//simulations_output/simulated_repall.vcf.WCAN1.ES17.FR.chr1.1PerCentMaf.3PCLRoutput", sep = "\t", header = T)
Set2.3_xpclr_chrall<- read.table(Set2.3_xpclr_chrall_new.txt", sep = "\t", header=T)

sim_WCAN1_ES17_FR$X3PCLR.Anc.std <- (sim_WCAN1_ES17_FR$X3PCLR.Anc-mean(sim_WCAN1_ES17_FR$X3PCLR.Anc))/sd(sim_WCAN1_ES17_FR$X3PCLR.Anc)
sim_WCAN1_ES17_FR$X3PCLR.A.std <- (sim_WCAN1_ES17_FR$X3PCLR.A-mean(sim_WCAN1_ES17_FR$X3PCLR.A))/sd(sim_WCAN1_ES17_FR$X3PCLR.A)
sim_WCAN1_ES17_FR$X3PCLR.B.std <- (sim_WCAN1_ES17_FR$X3PCLR.B-mean(sim_WCAN1_ES17_FR$X3PCLR.B))/sd(sim_WCAN1_ES17_FR$X3PCLR.B)


#WCAN1

hist(sim_WCAN1_ES17_FR$X3PCLR.A)

sort(Set2.3_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.05*length(Set2.3_xpclr_chrall$X3PCLR.A)]
sort(Set2.3_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.01*length(Set2.3_xpclr_chrall$X3PCLR.A)]
sort(Set2.3_xpclr_chrall$X3PCLR.A, decreasing=TRUE)[0.005*length(Set2.3_xpclr_chrall$X3PCLR.A)]

sort(Set2.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.05*length(Set2.3_xpclr_chrall$X3PCLR.A.std)]
sort(Set2.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.01*length(Set2.3_xpclr_chrall$X3PCLR.A.std)]
sort(Set2.3_xpclr_chrall$X3PCLR.A.std, decreasing=TRUE)[0.005*length(Set2.3_xpclr_chrall$X3PCLR.A.std)]

sort(sim_WCAN1_ES17_FR$X3PCLR.A, decreasing = T)[0.05*length(sim_wcan1_wcan2_CECAN$X3PCLR.A)]
sort(sim_WCAN1_ES17_FR$X3PCLR.A, decreasing = T)[0.01*length(sim_WCAN1_ES17_FR$X3PCLR.A)]
sort(sim_WCAN1_ES17_FR$X3PCLR.A, decreasing = T)[0.005*length(sim_WCAN1_ES17_FR$X3PCLR.A)]

sort(sim_WCAN1_ES17_FR$X3PCLR.A.std, decreasing = T)[0.05*length(sim_WCAN1_ES17_FR$X3PCLR.A.std)]
sort(sim_WCAN1_ES17_FR$X3PCLR.A.std, decreasing = T)[0.01*length(sim_WCAN1_ES17_FR$X3PCLR.A.std)]
sort(sim_WCAN1_ES17_FR$X3PCLR.A.std, decreasing = T)[0.005*length(sim_WCAN1_ES17_FR$X3PCLR.A.std)]

#ES17
hist(sim_WCAN1_ES17_FR$X3PCLR.B)

sort(Set2.3_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.05*length(Set2.3_xpclr_chrall$X3PCLR.B)]
sort(Set2.3_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.01*length(Set2.3_xpclr_chrall$X3PCLR.B)]
sort(Set2.3_xpclr_chrall$X3PCLR.B, decreasing=TRUE)[0.005*length(Set2.3_xpclr_chrall$X3PCLR.B)]

sort(Set2.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.05*length(Set2.3_xpclr_chrall$X3PCLR.B.std)]
sort(Set2.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.01*length(Set2.3_xpclr_chrall$X3PCLR.B.std)]
sort(Set2.3_xpclr_chrall$X3PCLR.B.std, decreasing=TRUE)[0.005*length(Set2.3_xpclr_chrall$X3PCLR.B.std)]

sort(sim_WCAN1_ES17_FR$X3PCLR.B, decreasing = T)[0.05*length(sim_WCAN1_ES17_FR$X3PCLR.B)]
sort(sim_WCAN1_ES17_FR$X3PCLR.B, decreasing = T)[0.01*length(sim_WCAN1_ES17_FR$X3PCLR.B)]
sort(sim_WCAN1_ES17_FR$X3PCLR.B, decreasing = T)[0.005*length(sim_WCAN1_ES17_FR$X3PCLR.B)]

sort(sim_WCAN1_ES17_FR$X3PCLR.B.std, decreasing = T)[0.05*length(sim_WCAN1_ES17_FR$X3PCLR.B.std)]
sort(sim_WCAN1_ES17_FR$X3PCLR.B.std, decreasing = T)[0.01*length(sim_WCAN1_ES17_FR$X3PCLR.B.std)]
sort(sim_WCAN1_ES17_FR$X3PCLR.B.std, decreasing = T)[0.005*length(sim_WCAN1_ES17_FR$X3PCLR.B.std)]

#ANC
hist(sim_WCAN1_ES17_FR$X3PCLR.Anc)

sort(Set2.3_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.05*length(Set2.3_xpclr_chrall$X3PCLR.Anc)]
sort(Set2.3_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.01*length(Set2.3_xpclr_chrall$X3PCLR.Anc)]
sort(Set2.3_xpclr_chrall$X3PCLR.Anc, decreasing=TRUE)[0.005*length(Set2.3_xpclr_chrall$X3PCLR.Anc)]

sort(Set2.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.05*length(Set2.3_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set2.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.01*length(Set2.3_xpclr_chrall$X3PCLR.Anc.std)]
sort(Set2.3_xpclr_chrall$X3PCLR.Anc.std, decreasing=TRUE)[0.005*length(Set2.3_xpclr_chrall$X3PCLR.Anc.std)]

sort(sim_WCAN1_ES17_FR$X3PCLR.Anc, decreasing = T)[0.05*length(sim_WCAN1_ES17_FR$X3PCLR.Anc)]
sort(sim_WCAN1_ES17_FR$X3PCLR.Anc, decreasing = T)[0.01*length(sim_WCAN1_ES17_FR$X3PCLR.Anc)]
sort(sim_WCAN1_ES17_FR$X3PCLR.Anc, decreasing = T)[0.005*length(sim_WCAN1_ES17_FR$X3PCLR.Anc)]

sort(sim_WCAN1_ES17_FR$X3PCLR.Anc.std, decreasing = T)[0.05*length(sim_WCAN1_ES17_FR$X3PCLR.Anc.std)]
sort(sim_WCAN1_ES17_FR$X3PCLR.Anc.std, decreasing = T)[0.01*length(sim_WCAN1_ES17_FR$X3PCLR.Anc.std)]
sort(sim_WCAN1_ES17_FR$X3PCLR.Anc.std, decreasing = T)[0.005*length(sim_WCAN1_ES17_FR$X3PCLR.Anc.std)]




###########QQ plots


###Set4

x <- sim_wcan1_wcan2_CECAN$X3PCLR.A
y <- Set3.1_xpclr_chrall$X3PCLR.B

x <- sim_wcan1_wcan2_CECAN$X3PCLR.B
y <- Set3.1_xpclr_chrall$X3PCLR.A

x <- sim_wcan1_wcan2_CECAN$X3PCLR.Anc
y <- Set3.1_xpclr_chrall$X3PCLR.Anc

x <- x[!is.na(x)]
y <- y[!is.na(y)]

x_std <- scale(x)[,1]
y_std <- scale(y)[,1]

n <- min(length(x_std), length(y_std))
probs <- seq(0, 1, length.out = n)

df_qq <- data.frame(
  qx = quantile(x_std, probs),
  qy = quantile(y_std, probs)
)

# ---- upper tail probabilities ----
tail_probs <- c(0.95, 0.99, 0.995)

df_tail <- data.frame(
  label = c("Top 5%", "Top 1%", "Top 0.5%"),
  qx = quantile(x_std, tail_probs),
  qy = quantile(y_std, tail_probs)
)

#For each set 3 plots (one for each branch)

p_s3.1_WCAN1<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(WCAN1)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s3.1_WCAN2<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +   
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(WCAN2)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s3.1_ANC<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="A(WCAN1, WCAN2)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))




#Set4

x <- sim_WCAN2_ES17_CECAN1$X3PCLR.A
y <- Set3.2_xpclr_chrall$X3PCLR.A

x <- sim_WCAN2_ES17_CECAN1$X3PCLR.B
y <- Set3.2_xpclr_chrall$X3PCLR.B

x <- sim_WCAN2_ES17_CECAN1$X3PCLR.Anc
y <- Set3.2_xpclr_chrall$X3PCLR.Anc

x <- x[!is.na(x)]
y <- y[!is.na(y)]

x_std <- scale(x)[,1]
y_std <- scale(y)[,1]

n <- min(length(x_std), length(y_std))
probs <- seq(0, 1, length.out = n)

df_qq <- data.frame(
  qx = quantile(x_std, probs),
  qy = quantile(y_std, probs)
)

# ---- upper tail probabilities ----
tail_probs <- c(0.95, 0.99, 0.995)

df_tail <- data.frame(
  label = c("Top 5%", "Top 1%", "Top 0.5%"),
  qx = quantile(x_std, tail_probs),
  qy = quantile(y_std, tail_probs)
)


p_s3.2_WCAN2<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(WCAN2)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s3.2_ES17 <- ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) + 
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(ES17)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s3.2_ANC<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="A(WCAN2, ES17)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))


#Set6

x <- sim_WCAN1_ES17_CECAN1$X3PCLR.A
y <- Set3.3_xpclr_chrall$X3PCLR.A

x <- sim_WCAN1_ES17_CECAN1$X3PCLR.B
y <- Set3.3_xpclr_chrall$X3PCLR.B

x <- sim_WCAN1_ES17_CECAN1$X3PCLR.Anc
y <- Set3.3_xpclr_chrall$X3PCLR.Anc

x <- x[!is.na(x)]
y <- y[!is.na(y)]

x_std <- scale(x)[,1]
y_std <- scale(y)[,1]

n <- min(length(x_std), length(y_std))
probs <- seq(0, 1, length.out = n)

df_qq <- data.frame(
  qx = quantile(x_std, probs),
  qy = quantile(y_std, probs)
)

# ---- upper tail probabilities ----
tail_probs <- c(0.95, 0.99, 0.995)

df_tail <- data.frame(
  label = c("Top 5%", "Top 1%", "Top 0.5%"),
  qx = quantile(x_std, tail_probs),
  qy = quantile(y_std, tail_probs)
)


p_s3.3_WCAN1<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(WCAN1)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s3.3_ES17 <- ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(ES17)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s3.3_ANC<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="A(WCAN1, ES17)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))



#Set1

x <- sim_wcan1_wcan2_fr$X3PCLR.A
y <- Set2.1_xpclr_chrall$X3PCLR.B

x <- sim_wcan1_wcan2_fr$X3PCLR.B
y <- Set2.1_xpclr_chrall$X3PCLR.A

x <- sim_wcan1_wcan2_fr$X3PCLR.Anc
y <- Set2.1_xpclr_chrall$X3PCLR.Anc

x <- x[!is.na(x)]
y <- y[!is.na(y)]

x_std <- scale(x)[,1]
y_std <- scale(y)[,1]

n <- min(length(x_std), length(y_std))
probs <- seq(0, 1, length.out = n)

df_qq <- data.frame(
  qx = quantile(x_std, probs),
  qy = quantile(y_std, probs)
)

# ---- upper tail probabilities ----
tail_probs <- c(0.95, 0.99, 0.995)

df_tail <- data.frame(
  label = c("Top 5%", "Top 1%", "Top 0.5%"),
  qx = quantile(x_std, tail_probs),
  qy = quantile(y_std, tail_probs)
)


p_s2.1_WCAN1<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(WCAN1)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s2.1_WCAN2<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5)  +
  scale_y_continuous(breaks=c(0,5,10))+
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(WCAN2)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s2.1_ANC<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) + 
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="A(WCAN1, WCAN2)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))




#Set2

x <- sim_WCAN2_ES17_FR$X3PCLR.A
y <- Set2.2_xpclr_chrall$X3PCLR.A

x <- sim_WCAN2_ES17_FR$X3PCLR.B
y <- Set2.2_xpclr_chrall$X3PCLR.B

x <- sim_WCAN2_ES17_FR$X3PCLR.Anc
y <- Set2.2_xpclr_chrall$X3PCLR.Anc

x <- x[!is.na(x)]
y <- y[!is.na(y)]

x_std <- scale(x)[,1]
y_std <- scale(y)[,1]

n <- min(length(x_std), length(y_std))
probs <- seq(0, 1, length.out = n)

df_qq <- data.frame(
  qx = quantile(x_std, probs),
  qy = quantile(y_std, probs)
)

# ---- upper tail probabilities ----
tail_probs <- c(0.95, 0.99, 0.995)

df_tail <- data.frame(
  label = c("Top 5%", "Top 1%", "Top 0.5%"),
  qx = quantile(x_std, tail_probs),
  qy = quantile(y_std, tail_probs)
)


p_s2.2_WCAN2<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(WCAN2)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s2.2_ES17<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(ES17)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s2.2_ANC<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="A(WCAN2, ES17)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))


#Set3

x <- sim_WCAN1_ES17_FR$X3PCLR.A
y <- Set2.3_xpclr_chrall$X3PCLR.A

x <- sim_WCAN1_ES17_FR$X3PCLR.B
y <- Set2.3_xpclr_chrall$X3PCLR.B

x <- sim_WCAN1_ES17_FR$X3PCLR.Anc
y <- Set2.3_xpclr_chrall$X3PCLR.Anc


x <- x[!is.na(x)]
y <- y[!is.na(y)]

x_std <- scale(x)[,1]
y_std <- scale(y)[,1]

n <- min(length(x_std), length(y_std))
probs <- seq(0, 1, length.out = n)

df_qq <- data.frame(
  qx = quantile(x_std, probs),
  qy = quantile(y_std, probs)
)

# ---- upper tail probabilities ----
tail_probs <- c(0.95, 0.99, 0.995)

df_tail <- data.frame(
  label = c("Top 5%", "Top 1%", "Top 0.5%"),
  qx = quantile(x_std, tail_probs),
  qy = quantile(y_std, tail_probs)
)

library(ggplot2)


#
p_s2.3_WCAN1<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(WCAN1)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s2.3_ES17<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) +
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="T(ES17)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))

p_s2.3_ANC<-ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail[1,], aes(qx, qy), color = "blue", size = 5, shape=15) +
  geom_point(data = df_tail[2,], aes(qx, qy), color = "green4", size = 5, shape=17) + 
  geom_point(data = df_tail[3,], aes(qx, qy), color = "orange2", size = 5, shape=19) +
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="A(WCAN1, ES17)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5))


#legend
df_tail$label <-as.factor(df_tail$label)
ggplot(df_qq, aes(qx, qy)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1) +
  geom_point(alpha = 0.5, size=2.5) +
  geom_point(data = df_tail, aes(x = qx,y =  qy, color = label, shape=label), size = 12) +
  scale_color_manual(values = c("orange2", "green4", "blue"))+
  scale_shape_manual(values = c(19, 17, 15))+
  # geom_text(data = df_tail,
  #           aes(qx, qy, label = label),
  #           hjust = -0.1,
  #           vjust = 2, size=8) +
  labs(x = "Simulated CLR-score ",
       y = "Observed CLR-score",
       title="A(WCAN1, ES17)")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(size=30),
        axis.title = element_blank(),
        plot.title = element_text(size = 30, hjust=0.5), 
        legend.key = element_rect(fill=NA),
        legend.key.size = unit(2,"cm"),
        legend.text = element_text(size=30)
        )



#Plot all together
library(cowplot)
plot_grid(plotlist = list(p_s2.1_WCAN2,p_s2.1_WCAN1, p_s2.1_ANC,
                          p_s2.2_WCAN2,p_s2.2_ES17, p_s2.2_ANC,
                          p_s2.3_WCAN1,p_s2.3_ES17, p_s2.3_ANC,
                          p_s3.1_WCAN2,p_s3.1_WCAN1, p_s3.1_ANC,
                          p_s3.2_WCAN2,p_s3.2_ES17, p_s3.2_ANC,
                          p_s3.3_WCAN1,p_s3.3_ES17, p_s3.3_ANC), ncol = 3, nrow=6)
