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
#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

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

## mafft per orthogroup
while read orthogroup; do
    echo ${orthogroup}
    mafft --maxiterate 1000 --globalpair --thread ${SLURM_NTASKS} ${SCO_DIR}/${orthogroup} > \
        ${basename}_msas/${orthogroup}.msa
done < sco_mapfile

module purge

# concatenate msas
/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/home/users/pierrj/git/python/concatenate_msas.py ${basename}_msas ${ACCESSIONS} ${basename}_msas_combined.fasta

module load trimal

# trim msa
trimal -gt 0.8 -in ${basename}_msas_combined.fasta -out ${basename}_msas_combined.trimmed.fasta