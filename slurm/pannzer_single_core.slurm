#!/bin/bash
#SBATCH --job-name=pannzer_single_core
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/scratch/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/scratch/users/pierrj/slurm_stderr/slurm-%j.out

## launch with
# proteome=something
# PROTEOMES_PATH=something

#  sbatch --job-name=$proteome.pannzer --mem=4000M -n 1 --export=ALL,proteome=$proteome,PROTEOMES_PATH=$PROTEOMES_PATH \
#     -p savio4_htc --qos=minium_htc4_normal --account=co_minium \
#     /global/home/users/pierrj/git/slurm/pannzer_single_core.slurm

source activate /global/scratch/users/pierrj/conda_envs/pannzer

cd $SLURM_SUBMIT_DIR

protein_count=$(grep gene GO.${proteome}.out | awk '{print $1}' | sort | uniq | wc -l)
b=7000

if (( protein_count < b )); then
    mv GO.${proteome}.out GO.${proteome}.out.old 
    /global/scratch/users/pierrj/conda_envs/pannzer/bin/python \
        /global/scratch/users/pierrj/pannzer2/SANSPANZ.3/runsanspanz.py \
        -R -o ",DE.${proteome}.out,GO.${proteome}.out,anno.${proteome}.out" -s "Pyricularia oryzae" < ${PROTEOMES_PATH}/${proteome}
fi