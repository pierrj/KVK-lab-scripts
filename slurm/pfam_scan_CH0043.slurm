#!/bin/bash
#SBATCH --job-name=pfam_scan_CH0043
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes/pfam_scan

source activate /global/scratch/users/pierrj/conda_envs/pfam_scan

OUTPUTNAME=CH0043
THREADS=${SLURM_NTASKS}

pfam_scan.pl -cpu ${THREADS} -outfile ${OUTPUT_NAME}.pfamscan.out -dir ${PFAM_DIR} -fasta ${OUTPUT_NAME}.CDS.fasta

/global/scratch/users/pierrj/scripts/plant_rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl -p ${OUTPUT_NAME}.pfamscan.out -e 0.001 -o ${OUTPUT_NAME}.pfamscan.kparse.out
