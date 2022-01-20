#!/bin/bash
#SBATCH --job-name=effector_predictions_gladieux_proteomes
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/moryzae_pav/effector_annotations

PROTEOMES_PATH=/global/scratch/users/pierrj/moryzae_pav/gladieux_et_al_2020_data/gladieux_et_al_2021_annotations
MAPFILE=gladieux_et_al_2021_assemblies_mapfile

module load seqtk

while read genome; do
    /global/home/users/pierrj/git/bash/effector_predictor.sh -i ${PROTEOMES_PATH}/${genome}/${genome}_protein.fasta -o ${genome}
do < $MAPFILE