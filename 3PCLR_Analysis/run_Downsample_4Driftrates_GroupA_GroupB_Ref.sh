#!/bin/bash

java -jar DownsampleSFS_final_bugfixed.jar\
	--vcf Cantabrian_GroupA_GroupB_Ref_10Permissing_polarized_commonsites_CanParaOut_FranParaOut_noGenes2kb.b.vcf.gz\
	--cutoff 0.1\
	--popfile $popfile\
	--listofpop $GroupA,$GroupB,$Ref\
	--haploidize
