#!/bin/bash
#SBATCH --job-name=public_github_test
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_caller_test_post_comments/RC_1B

export ECC_CALLER_PYTHON_SCRIPTS=/global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_caller_test_post_comments/ecc_caller/python_scripts

export PATH=$PATH:/global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_caller_test_post_comments/ecc_caller

mapfile="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_no_organelles.contignames"
genome_bwa="/global/scratch/users/pierrj/references/ORSA_IRGSP-1.0_bwa"

# generate_bam_file.sh -g $genome_bwa -1 RC_1B_R1.fastq -2 RC_1B_R2.fastq -s RC_1B -t ${SLURM_NTASKS} -m $mapfile

# call_ecc_regions.sh -m $mapfile -s RC_1B -t ${SLURM_NTASKS} -b filtered.sorted.RC_1B.bam

assign_confidence.sh -m $mapfile -s RC_1B -t ${SLURM_NTASKS} -b filtered.sorted.RC_1B.bam -r RC_1B.confirmedsplitreads.bed
