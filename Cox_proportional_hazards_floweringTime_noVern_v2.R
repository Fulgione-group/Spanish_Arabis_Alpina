
library(survival)
library(survminer)
library(ggplot2)

# cox proportional hazards-model on the flowering time data without vernalization



data <- read.csv("input_data_flowering_haplotype.csv", header = T)


# the experiment ended at 250 DAS
  # flowering / non-flowering = 1/0
# haplotype 1 = 0
# haplotype 2 = 2



# haplogroup as factor
data$haplotype <- factor(data$haplotype)


res.cox <- coxph(Surv(DAS, flowering) ~ haplotype, data = data)


# model with only haplogroup as predictor 
 
summary(res.cox)


#coxph(formula = Surv(DAS, flowering) ~ haplogroup, data = data)

#n= 479, number of events= 216 

#coef exp(coef) se(coef)      z Pr(>|z|)    
#haplotype2 -3.15323   0.04271  0.16831 -18.73   <2e-16 ***
 # ---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#exp(coef) exp(-coef) lower .95 upper .95
#haplogroup2   0.04271      23.41   0.03071   0.05941

#Concordance= 0.799  (se = 0.009 )
#Likelihood ratio test= 397.4  on 1 df,   p=<2e-16
#Wald test            = 351  on 1 df,   p=<2e-16
#Score (logrank) test = 567.3  on 1 df,   p=<2e-16





# plot the baseline survival function


ggsurvplot(survfit(res.cox), palette = "#1F77B4", data = data, 
           ggtheme = theme_minimal())


#ggsave("baseline_survival_function.png", width = 12, height = 6)



sex_df <- with(data,
               data.frame(haplotype = c("0", "2")
               )
)
sex_df



fit <- survfit(res.cox, newdata = sex_df)


ggsurvplot(fit, conf.int = TRUE, legend.title=" ", legend.labs=c("Haplotype 1", "Haplotype 2"), 
           data = data, 
           ylab=c("proportion not flowering", cex=2), 
           xlab=c("days from germination", cex=2), 
           palette=c("#1F77B4", "#FFC300"),
           ggtheme = theme_minimal() +
             theme(legend.text = element_text(size = 14),
                   axis.text = element_text(size = 18), 
                   axis.title = element_text(size = 18)))



#ggsave("haplotypes_model_2.pdf", width = 14, height = 8)





# include region as a covariate

head(data)

res.cox <- coxph(Surv(DAS, flowering) ~ haplotype + region, data = data)

summary(res.cox)


#n= 479, number of events= 216 

#coef exp(coef) se(coef)      z Pr(>|z|)    
#haplotype2 -1.26092   0.28339  0.17017  -7.41 1.27e-13 ***
 # regionIBE   -3.54964   0.02873  0.32560 -10.90  < 2e-16 ***
  #---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#exp(coef) exp(-coef) lower .95 upper .95
#haplogroup2   0.28339      3.529   0.20302    0.3956
#regionIBE     0.02873     34.801   0.01518    0.0544

#Concordance= 0.879  (se = 0.008 )
#Likelihood ratio test= 558.4  on 2 df,   p=<2e-16
#Wald test            = 259.9  on 2 df,   p=<2e-16
#Score (logrank) test = 657.2  on 2 df,   p=<2e-16



