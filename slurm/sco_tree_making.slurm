#!/bin/bash
#SBATCH --job-name=sco_tree_making
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/moryzae_pav/tree_making

ls Single_Copy_Orthologue_Sequences_jan_04_2021 > sco_mapfile

if [ -d "Single_Copy_Orthologue_Sequences_jan_04_2021_msas" ]; then
    rm -r Single_Copy_Orthologue_Sequences_jan_04_2021_msas
fi

while read orthogroup; do
    mafft --maxiterate 1000 --globalpair --thread 24 Single_Copy_Orthologue_Sequences_jan_04_2021/${orthogroup} > \
        Single_Copy_Orthologue_Sequences_jan_04_2021_msas/${orthogroup}.msa
done < sco_mapfile