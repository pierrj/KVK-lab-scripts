#!/bin/bash
#SBATCH --job-name=local_copy_time_comparison
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
cd /local/job${SLURM_JOBID}/

cp /global/scratch/users/pierrj/references/ORSA_w_organelles_bwa.amb /local/job${SLURM_JOBID}/

cp /global/scratch/users/pierrj/references/ORSA_w_organelles_bwa.ann /local/job${SLURM_JOBID}/

cp /global/scratch/users/pierrj/references/ORSA_w_organelles_bwa.bwt /local/job${SLURM_JOBID}/

cp /global/scratch/users/pierrj/references/ORSA_w_organelles_bwa.pac /local/job${SLURM_JOBID}/

cp /global/scratch/users/pierrj/references/ORSA_w_organelles_bwa.sa /local/job${SLURM_JOBID}/

cp /global/scratch/users/pierrj/references/ORSA_w_organelles.fasta /local/job${SLURM_JOBID}/

cp /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A_R1.fastq /local/job${SLURM_JOBID}/

cp /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A_R2.fastq /local/job${SLURM_JOBID}/

date +"%T"

bwa mem -t 24 ORSA_w_organelles_bwa RC_1A_R1.fastq RC_1A_R2.fastq -o tmp.bwamem_local_test

date +"%T"

rm tmp.bwamem_local_test


cd /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A

date +"%T"

bwa mem -t 24 /global/scratch/users/pierrj/references/ORSA_w_organelles_bwa RC_1A_R1.fastq RC_1A_R2.fastq -o tmp.bwamem_local_test

date +"%T"

rm tmp.bwamem_local_test