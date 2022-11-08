#!/bin/bash
#SBATCH --job-name=pannzer_per_proteome
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --cpus=1
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

DATA_OUTPUT_DIR=go

cd $PROJ_DIR

if [ -d $DATA_OUTPUT_DIR ]; then
    rm -r $DATA_OUTPUT_DIR
fi

mkdir $DATA_OUTPUT_DIR

cd $DATA_OUTPUT_DIR

source activate /global/scratch/users/pierrj/conda_envs/pannzer

/global/scratch/users/pierrj/conda_envs/pannzer/bin/python \
    /global/scratch/users/pierrj/pannzer2/SANSPANZ.3/runsanspanz.py \
    -R -o ",DE.out,GO.out,anno.out" -s "Pyricularia oryzae" < $PROTEIN_FASTA