#!/bin/bash
#SBATCH --job-name=effector_predictions
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/predicted_effectors

PROTEOMES_PATH=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/all_proteomes/
MAPFILE=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/genomes_mapfile

module load seqtk

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

while read proteome; do
    echo /global/home/users/pierrj/git/bash/effector_predictor.sh -i ${PROTEOMES_PATH}/${proteome}_fungap_out_prot.faa -o ${proteome}.effectors_pred >> jobqueue
done < ${MAPFILE}

parallel -j ${SLURM_NTASKS} < jobqueue