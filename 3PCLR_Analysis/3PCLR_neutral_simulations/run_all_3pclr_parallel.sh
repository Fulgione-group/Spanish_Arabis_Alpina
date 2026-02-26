#!/bin/bash
set -euo pipefail

source /netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/anaconda3/etc/profile.d/conda.sh
conda activate msprime-env

cd Sim_3PCLR

MAX_JOBS=100


run_triplet () {
	local P1=$1 P2=$2 O=$3 DRIFT=$4 VCF=$5

	bash run_CreatInput_3PCLR_SimDem_var.sh "$P1" "$P2" "$O" "$VCF"

	bash filterout_maf0.01.sh \
		"${VCF}.${P1}.${P2}.${O}.chr1..3PCLRinput" \
		"${VCF}.${P1}.${P2}.${O}.chr1.1PerCentMaf.3PCLRinput"

	bash run_3pclr_var.sh \
		"${VCF}.${P1}.${P2}.${O}.chr1.1PerCentMaf.3PCLRinput" \
		"${VCF}.${P1}.${P2}.${O}.chr1.1PerCentMaf.3PCLRoutput" \
		"$DRIFT"
}



run_rep () {
    	i=$1
    	echo ">>> starting replicate ${i}"

	python demography_3pclr_var.py \
  	--out-vcf "simulated_rep${i}.vcf" --rep "${i}"

	VCF="simulated_rep${i}.vcf"
	# Run 3PCLR with fixed drift rates from empirical analysis
	# FR outgroup (3 jobs)
	run_triplet WCAN1 WCAN2 FR   "0.165787065311824,0.176538266986559,3.35270613270779" "$VCF" &
	run_triplet WCAN1 ES17  FR   "0.637915215732256,0.309793863105167,2.84623750463338" "$VCF" &
	run_triplet WCAN2 ES17  FR   "0.653330723082598,0.31866019466855,2.77826372683411" "$VCF" &

	# CECAN outgroup (3 jobs)
	run_triplet WCAN1 WCAN2 CECAN  "0.0894536733072693,0.0823830737100736,0.414595628164814" "$VCF" &
	run_triplet WCAN1 ES17  CECAN1 "0.218110931996827,0.249701987575031,0.41869017875587"   "$VCF" &
	run_triplet WCAN2 ES17  CECAN1 "0.223249501776188,0.264449526761651,0.408417855269791"  "$VCF" &

	wait   # wait for the 6 pipelines of this replicate

    	echo "<<< finished replicate ${i}"
}

for i in $(seq 1 80); do
    # throttle to MAX_JOBS
    while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]; do
        sleep 5
    done

    run_rep "${i}" > log_rep${i}.out 2>&1 &
done

wait
echo "ALL REPLICATES FINISHED"
