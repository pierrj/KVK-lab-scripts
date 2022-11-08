#!/bin/bash
#SBATCH --job-name=effector_predictor
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --cpus=1
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

DATA_OUTPUT_DIR=effector_predictor

cd $PROJ_DIR

if [ -d $DATA_OUTPUT_DIR ]; then
    rm -r $DATA_OUTPUT_DIR
fi

mkdir $DATA_OUTPUT_DIR

cd $DATA_OUTPUT_DIR

/global/home/users/pierrj/git/bash/effector_predictor.sh -i $PROTEIN_FASTA -o output