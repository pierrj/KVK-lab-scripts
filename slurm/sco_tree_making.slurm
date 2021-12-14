#!/bin/bash
#SBATCH --job-name=sco_tree_making
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/moryzae_pav/tree_making

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

# cp -r Single_Copy_Orthologue_Sequences Single_Copy_Orthologue_Sequences_dec_08_21

SCO_DIR=Single_Copy_Orthologue_Sequences_dec_08_21

# ls ${SCO_DIR} > sco_mapfile

# if [ -d "${SCO_DIR}_msas" ]; then
#     rm -r ${SCO_DIR}_msas
# fi

# mkdir ${SCO_DIR}_msas

# while read orthogroup; do
#     mafft --maxiterate 1000 --globalpair --thread ${SLURM_NTASKS} ${SCO_DIR}/${orthogroup} > \
#         ${SCO_DIR}_msas/${orthogroup}.msa
# done < sco_mapfile

# module remove python

# python /global/home/users/pierrj/git/python/concatenate_msas.py ${SCO_DIR}_msas proteomes_combined_for_orthofinder_accessions_re ${SCO_DIR}_msas_combined.fasta

module load trimal

trimal -gt 0.8 -in ${SCO_DIR}_msas_combined.fasta -out ${SCO_DIR}_msas_combined.trimmed.fasta

fasttree -gamma -out ${SCO_DIR}_msas_combined.trimmed.fasttree ${SCO_DIR}_msas_combined.trimmed.fasta