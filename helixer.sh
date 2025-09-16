#!/bin/bash -l
#SBATCH -o ./helixer_%j.out
#SBATCH -e ./helixer_%j.err
#SBATCH -D ./
#SBATCH -J helixer
#SBATCH --constraint="gpu"
#SBATCH --gres=gpu:a100:1       
#SBATCH --cpus-per-task=18
#SBATCH --mem=125000
#SBATCH --ntasks=1
#SBATCH --time=24:00:00


ml purge
ml load intel/21.2.0 impi/2021.2 cuda/11.2
ml singularity

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

OUTDIR=
cd $OUTDIR
pwd -P
FASTA=Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz

SIF=helixer-docker_helixer_v0.3.0a0_cuda_11.2.0-cudnn8.sif

singularity run --nv --bind $OUTDIR $SIF Helixer.py --fasta-path $FASTA --lineage land_plant --gff-output-path $OUTDIR/Arabis_alpina.MPIPZ.version_5.1_helixer.gff3
