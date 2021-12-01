#!/bin/bash
#SBATCH --job-name=sco_tree_making
#SBATCH --partition=savio_bigmem
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

# ls Single_Copy_Orthologue_Sequences_nov_01_2021 > sco_mapfile

# if [ -d "Single_Copy_Orthologue_Sequences_nov_01_2021_msas" ]; then
#     rm -r Single_Copy_Orthologue_Sequences_nov_01_2021_msas
# fi

# mkdir Single_Copy_Orthologue_Sequences_nov_01_2021_msas

# while read orthogroup; do
#     mafft --maxiterate 1000 --globalpair --thread 24 Single_Copy_Orthologue_Sequences_nov_01_2021/${orthogroup} > \
#         Single_Copy_Orthologue_Sequences_nov_01_2021_msas/${orthogroup}.msa
# done < sco_mapfile

# module remove python

# python /global/home/users/pierrj/git/python/concatenate_msas.py Single_Copy_Orthologue_Sequences_nov_01_2021_msas proteomes_combined_for_orthofinder_accessions SCO_nov_01_2021_combined_msas.fasta

# module load trimal

# trimal -gt 0.8 -in SCO_nov_01_2021_combined_msas.fasta -out SCO_nov_01_2021_combined_msas.trimmed.fasta

fasttree -gamma -out SCO_nov_01_2021_combined_msas.trimmed.fasttree SCO_nov_01_2021_combined_msas.trimmed.fasta