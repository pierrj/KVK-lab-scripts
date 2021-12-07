#!/bin/bash
#SBATCH --job-name=orthogroup_consensus_sequences
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/moryzae_pav/orthogroup_cons


source activate /global/scratch/users/pierrj/conda_envs/orthofinder

module load trimal
module load emboss

ls -1 orthogroup_cds_out > orthogroup_mapfile

input_dir=orthogroup_protein_out
out_msas=orthogroup_protein_out_msas
out_cons=orthogroup_protein_out_cons

if [ -d "${out_msas}" ]; then
    rm -r ${out_msas}
fi

mkdir ${out_msas}

if [ -d "${out_cons}" ]; then
    rm -r ${out_cons}
fi

mkdir ${out_cons}

while read orthogroup; do
    mafft --maxiterate 1000 --globalpair --thread 24 ${input_dir}/${orthogroup} > \
        ${out_msas}/${orthogroup}.msa
    trimal -gt 0.8 -in ${out_msas}/${orthogroup}.msa -out ${out_msas}/${orthogroup}.trimmed.msa
    cons -sequence ${out_msas}/${orthogroup}.msa -out ${out_cons}/${orthogroup}.cons
done < sco_mapfile