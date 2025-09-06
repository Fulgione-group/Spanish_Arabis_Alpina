# infer multicopy regions in Cantabrian and French populations


## Cantabria

Prepare input from VCF with max 30% missing genotypes per site:  "run_PrepareParaMaskInput_fromVCF_Cantabria_v2.sh"
Run ParaMask EM step: "run_ParaMask_Cantabria_v2_2.7.1.sh"
Cluster seed snps per chromosome: "run_ParaMask_Cluster_Seeds_Cantabrian_final_chr*.sh"

## France

Prepare input from VCF with max 30% missing genotypes per site:  "run_PrepareParaMaskInput_French_final_v2.sh"
Run ParaMask EM step: "run_ParaMask_EM_French_clean_2.7.1.sh"
Cluster seed snps per chromosome: "run_ParaMask_Cluster_Seeds_French_chr*.sh"

