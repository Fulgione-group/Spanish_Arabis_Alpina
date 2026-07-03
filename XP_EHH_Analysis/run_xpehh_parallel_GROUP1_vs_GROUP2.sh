mkdir -p logs results

CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8)
CHR=${CHRS[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "$(date) Starting XP-EHH for ${CHR}"

selscan \
  --xpehh \
  --trunc-ok \
  --max-gap 500000 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --vcf GROUP1.${CHR}.vcf.gz \
  --vcf-ref GROUP2.${CHR}.vcf.gz \
  --map selscan.${CHR}.map \
  --out results/GROUP1_vs_GROUP2.${CHR}

echo "$(date) Finished XP-EHH for ${CHR}"
