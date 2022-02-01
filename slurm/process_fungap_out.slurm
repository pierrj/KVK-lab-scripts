#!/bin/bash
#SBATCH --job-name=process_fungap_out
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/full_run_orthofinder_1_31_2022

if [ ! -d "full_runfungap_proteomes_1_31_2022" ]; then
    mkdir full_runfungap_proteomes_1_31_2022
fi

while read genome; do
    cp /global/scratch/users/pierrj/fungap_runs/gladieux_all/${genome}/fungap_out/fungap_out/fungap_out_prot.faa \
        full_runfungap_proteomes_1_31_2022/${genome}_fungap_out_prot.faa
done < /global/scratch/users/pierrj/fungap_runs/gladieux_all/genomes_mapfile

python /global/home/users/pierrj/git/python/process_protein_sequences_for_orthofinder.py full_runfungap_proteomes_1_31_2022 \
    full_runfungap_proteomes_1_31_2022_corrected \
    gladieux_et_al_assembly_stats.txt \
    output_accession_names