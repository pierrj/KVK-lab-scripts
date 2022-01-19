#!/bin/bash
#SBATCH --job-name=plot_profiler_test
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/plot_profiler


#awk -F "\t|;" -v OFS='\t' '{ if ( $5 == "<DEL>") {print $1, $2, substr($14, 5, length($14)-5}}' ERR4832458.all.vcf

REGIONS_BED=ERR4832458.dels.bed
DENSITY_BED=moryzae_70-15_ref_with_mito.fasta.out.gff
GENOME_FILE=moryzae_70-15_ref_with_mito.fasta
WINDOWS=10
THREADS=${SLURM_NTASKS}
OUTPUT_NAME=test

source activate /global/scratch/users/pierrj/conda_envs/deeptools

/global/home/users/pierrj/git/bash/plot_profiler.sh -b $REGIONS_BED -d $DENSITY_BED -g $GENOME_FILE -w $WINDOWS -t $THREADS -o $OUTPUT_NAME