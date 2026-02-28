# Test genotype associations of FRL1 with flowering time (without vernalization) and NAC055 with aridity index

NAC055 is tested within Cantabria, and FRL1 is tested in Cantabria and France together.  
Structure covariates are included in both analyses.

- Binomial GLM using the base R `glm()` function, with NAC055 insertion allele dosage as predictor and response:  
  `NAC055_AI_corr.R`

- Linear mixed model using the `lme4` package (v1.1-28), modeling replicates of genotypes as random effect:  
  `Flowering_FRL1_v2.R`

- Cox proportional hazard (of flowering) model with and without region as covariate:
  `Cox_proportional_hazards_floweringTime_noVern_v2.R`
