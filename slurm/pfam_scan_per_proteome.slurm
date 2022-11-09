#!/bin/bash
#SBATCH --job-name=pfam_scan_per_proteome
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

DATA_OUTPUT_DIR=pfam_out

cd $PROJ_DIR

if [ -d $DATA_OUTPUT_DIR ]; then
    rm -r $DATA_OUTPUT_DIR
fi

mkdir $DATA_OUTPUT_DIR

cd $DATA_OUTPUT_DIR

source activate /global/scratch/users/pierrj/conda_envs/pfam_scan

pfam_scan.pl -cpu ${SLURM_NTASKS} -outfile pfamscan.out -dir ${PFAM_DIR} -fasta $PROTEIN_FASTA

/global/scratch/users/pierrj/parse_pfam_domains/plant_rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl -p pfamscan.out -e 0.001 -o pfamscan.kparse.out

awk -v OFS='\t' '{print $1, "TRUE"}' pfamscan.kparse.out  > pfamscan.table