#!/bin/bash
#SBATCH --job-name=transfer_guy11_info
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

DATA_OUTPUT_DIR=transfer_guy11

cd $PROJ_DIR

if [ -d $DATA_OUTPUT_DIR ]; then
    rm -r $DATA_OUTPUT_DIR
fi

mkdir $DATA_OUTPUT_DIR

cd $DATA_OUTPUT_DIR

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

module purge

export PERL5LIB=''

mkdir guy11_plus_input

cp $GUY11_PROTEOME guy11_plus_input/
cp $PROTEIN_FASTA guy11_plus_input/

orthofinder -S diamond_ultra_sens -f guy11_plus_input -M msa -A mafft -T fasttree -t ${SLURM_NTASKS} -a ${SLURM_NTASKS} -o orthofinder_out

orthofinder_date=$(date +%b%d)

cp orthofinder_out/Results_${orthofinder_date}/Orthogroups/Orthogroups.tsv .

python /global/home/users/pierrj/git/python/transfer_info_from_guy11.py Orthogroups.tsv $GUY11_DATA_DIR