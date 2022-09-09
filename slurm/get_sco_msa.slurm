#!/bin/bash
#SBATCH --job-name=sco_msas
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

PROJECT_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all
DATE=Sep08

cd $PROJECT_DIR

SCO_DIR=$PROJECT_DIR/orthofinder/Results_${DATE}/WorkingDirectory/OrthoFinder/Results_${DATE}/Single_Copy_Orthologue_Sequences
ACCESSIONS=$PROJECT_DIR/output_accession_names

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

/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/concatenate_msas.py ${basename}_msas ${ACCESSIONS} ${basename}_msas_combined.fasta

module load trimal

trimal -gt 0.8 -in ${basename}_msas_combined.fasta -out ${basename}_msas_combined.trimmed.fasta