#!/bin/bash
#SBATCH --job-name=orthofinder_add_70_15
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/adding_70_15_to_orthogrouping

# cp -r ../full_run_orthofinder_1_31_2022/orthofinder_1_31_2022/ .

# mkdir 70_15_proteome
# cd 70_15_proteome
# wget http://ftp.ensemblgenomes.org/pub/fungi/release-52/fasta/magnaporthe_oryzae/pep/Magnaporthe_oryzae.MG8.pep.all.fa.gz
# gunzip Magnaporthe_oryzae.MG8.pep.all.fa.gz
# cd ..

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

module purge

export PERL5LIB=''

orthofinder -oa -M msa -A mafft -T fasttree -t ${SLURM_NTASKS} -a 24 -b orthofinder_1_31_2022/Results_Feb03/WorkingDirectory/ -f 70_15_proteome