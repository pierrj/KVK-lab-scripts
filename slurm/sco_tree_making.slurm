#!/bin/bash
#SBATCH --job-name=sco_msas
#SBATCH --partition=savio_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/full_run_tree_1_31_2022

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

# cp -r Single_Copy_Orthologue_Sequences Single_Copy_Orthologue_Sequences_dec_08_21

SCO_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/full_run_orthofinder_1_31_2022/orthofinder_1_31_2022/Results_Feb03/WorkingDirectory/OrthoFinder/Results_Feb03/Single_Copy_Orthologue_Sequences
ACCESSIONS=/global/scratch/users/pierrj/PAV_SV/PAV/full_run_orthofinder_1_31_2022/output_accession_names

basename=$(basename ${SCO_DIR})

ls ${SCO_DIR} > sco_mapfile

if [ -d "${basename}_msas" ]; then
    rm -r ${basename}_msas
fi

mkdir ${basename}_msas

while read orthogroup; do
    echo ${orthogroup}
    mafft --maxiterate 1000 --globalpair --thread ${SLURM_NTASKS} ${SCO_DIR}/${orthogroup} > \
        ${basename}_msas/${orthogroup}.msa
done < sco_mapfile

module purge

python /global/home/users/pierrj/git/python/concatenate_msas.py ${basename}_msas ${ACCESSIONS} ${basename}_msas_combined.fasta

module load trimal

trimal -gt 0.8 -in ${basename}_msas_combined.fasta -out ${basename}_msas_combined.trimmed.fasta

# fasttree -gamma -out ${SCO_DIR}_msas_combined.trimmed.fasttree ${SCO_DIR}_msas_combined.trimmed.fasta