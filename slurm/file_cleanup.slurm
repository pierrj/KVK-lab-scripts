#!/bin/bash
#SBATCH --job-name=file_cleanup
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

rm -r /global/scratch/users/pierrj/foldseek/mg8_proteome_tests
rm -r /global/home/users/pierrj/foldseek/mg8_proteome_tests

cd /global/scratch/users/pierrj/foldseek/

mkdir mg8_proteome_tests

cd mg8_proteome_tests

mkdir structures

cd structures

cp /global/scratch/users/pierrj/foldseek/proteome-tax_id-242507-0_v3.tar .
cp /global/scratch/users/pierrj/foldseek/proteome-tax_id-242507-1_v3.tar .

tar -xvf proteome-tax_id-242507-0_v3.tar
tar -xvf proteome-tax_id-242507-1_v3.tar

rm *.json.gz

rm proteome-tax_id-242507-0_v3.tar
rm proteome-tax_id-242507-1_v3.tar