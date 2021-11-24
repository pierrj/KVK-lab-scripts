#!/bin/bash
#SBATCH --job-name=ecc_caller_git_test
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


cd /global/scratch/users/pierrj/eccDNA/ecc_caller_git_test


module purge

source activate /global/scratch/users/pierrj/ecc_caller_test_2

export PERL5LIB=/global/scratch/users/pierrj/ecc_caller_test_2/lib


export ECC_CALLER_PYTHON_SCRIPTS=/global/scratch/users/pierrj/eccDNA/ecc_caller_git_test/ecc_caller/python_scripts/
export ECC_CALLER_PICARD=/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/
export PATH="/global/scratch/users/pierrj/eccDNA/ecc_caller_git_test/ecc_caller/:$PATH"

# bwa index -p genome_bwa genome.fasta 

# generate_bam_file.sh -g genome_bwa -1 R1.fastq -2 R2.fastq -s output_name -t 24 -m mapfile

call_ecc_regions.sh -m mapfile -s output_name -t 24 -b uniq.filtered.sorted.output_name.bam -q multimapped.filtered.name_sorted.output_name.bam

assign_confidence.sh -m mapfile -s output_name -t 24 -b no_secondary.filtered.sorted.output_name.bam -r output_name.confirmedsplitreads.bed