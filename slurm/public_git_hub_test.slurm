#!/bin/bash
#SBATCH --job-name=public_git_hub_test
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_caller_public_version_tests

cd ecc_caller
git pull
cd ..

export ECC_CALLER_PYTHON_SCRIPTS=ecc_caller/python_scripts/
export ECC_CALLER_PICARD=/clusterfs/vector/home/groups/software/sl-7.x86_64/modules/picard/2.9.0/lib/

/ecc_caller/generate_bam_file.sh -g genome_bwa -1 R1.fastq -2 R2.fastq -s output_name -t n_threads -m mapfile

/ecc_caller/call_ecc_regions.sh -m mapfile -s output_name -t n_threads -b uniq.filtered.sorted.output_name.bam -q multimapped.filtered.name_sorted.output_name.bam

/ecc_caller/assign_confidence.sh -m mapfile -s output_name -t n_threads -b no_secondary.filtered.sorted.output_name.bam -r output_name.confirmedsplitreads.bed