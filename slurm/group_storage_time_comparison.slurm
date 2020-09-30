#!/bin/bash
#SBATCH --job-name=greoupst_time_comparison
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/group_storage_time_comparison

date +"%T"

cp /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A_R1.fastq .

cp /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A_R2.fastq .

date +"%T"

cd /global/home/groups/fc_kvkallow/group_storage_time_comparison

date +"%T"

cp /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A_R1.fastq .

cp /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/RC_1A/RC_1A_R1.fastq .

date +"%T"


cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/group_storage_time_comparison

date +"%T"

bwa mem -t 24 /global/scratch/users/pierrj/references/ORSA_w_organelles_bwa RC_1A_R1.fastq RC_1A_R2.fastq -o tmp.bwamem_local_test

date +"%T"

rm tmp.bwamem_local_test

cd /global/home/groups/fc_kvkallow/group_storage_time_comparison

date +"%T"

bwa mem -t 24 /global/scratch/users/pierrj/references/ORSA_w_organelles_bwa RC_1A_R1.fastq RC_1A_R2.fastq -o tmp.bwamem_local_test

date +"%T"

rm tmp.bwamem_local_test