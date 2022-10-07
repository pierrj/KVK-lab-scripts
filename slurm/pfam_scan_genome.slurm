#!/bin/bash
#SBATCH --job-name=pfam_scan_genome
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

## launch with
# genome=something
# proteome=something_else
# sbatch --export=ALL,genome=$genome,proteome=$proteome /global/home/users/pierrj/git/slurm/histone_marks_per_genome.slurm


PROJ_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/${genome}_pfam_scan
THREADS=$SLURM_NTASKS

if [ -d "$PROJ_DIR" ]; then
    rm -r "$PROJ_DIR"
fi

mkdir $PROJ_DIR

cd $PROJ_DIR

PROTEOMES_PATH=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/all_proteomes_corrected/
FASTA_FILE=${PROTEOMES_PATH}/${proteome}
PFAM_DIR=/global/scratch/users/pierrj/pfam/

source activate /global/scratch/users/pierrj/conda_envs/pfam_scan

pfam_scan.pl -cpu ${THREADS} -outfile ${genome}.pfamscan.out -dir ${PFAM_DIR} -fasta $FASTA

/global/scratch/users/pierrj/parse_pfam_domains/plant_rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl -p ${genome}.pfamscan.out -e 0.001 -o ${genome}.pfamscan.kparse.out