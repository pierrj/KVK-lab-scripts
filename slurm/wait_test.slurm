#!/bin/bash
#SBATCH --job-name=wait_test
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


cd /global/scratch/users/pierrj/sandbox/wait_test

date '+%d/%m/%Y %H:%M:%S' >> master

if [ -f slurm_ids ]; then
    rm slurm_ids
fi

for((i=0; i<5; i++)); do
    sbatch --parsable --job-name=$i.wait_test --export=print=$i /global/home/users/pierrj/git/slurm/wait_command.slurm >> slurm_ids
done

sbatch -W --dependency=afterok:$(echo slurm_ids | tr ' ' :) /global/home/users/pierrj/git/slurm/wait_command.slurm
wait

date '+%d/%m/%Y %H:%M:%S' >> master