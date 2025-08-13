#3PCLR analysis

java source files and bytecode are in separate directories  "java_src" and "java_bytecode"

##Prepare input

Create input from a vcf per chromosome using a genetic map, missingness cutoff of 10%, to subject populations (Group1, Group2) and one reference population (REF): "run_CreatInput_3PCLR_chrX_G1_G2_Ref_v3.sh"


##3PCLR

run 3PCLR analysis per chromosome with minimum distance of 25 bp , 100 snps per window, : "run_3pclr_G1_G2_REF_chrX.sh"
