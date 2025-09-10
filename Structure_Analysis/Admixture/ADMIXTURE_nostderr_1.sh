#!/bin/bash

BED="CantabrianParaOut10PercentMissingnessNoSingletonsUnrelated_pruned.bed"
for k in {1..20}
do
	admixture --cv $BED $k -j3 &
done
wait
