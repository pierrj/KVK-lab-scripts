#!/bin/bash
#SBATCH --job-name=test_ecc_finder
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/test_ecc_sv_finder_script

REFERENCE=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
QUERY=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/GCA_000002495.2_MG8_genomic.fna
THREADS=${SLURM_NTASKS}
OUTPUT_NAME=test
PERCENT_ZEROES_FILTER=1.0
OUTPUT_DIR=output_jpgs

source ~/.bashrc
conda activate mummer

module load gnuplot
module load imagemagick

/global/home/users/pierrj/git/bash/mummerplotter.sh -r ${REFERENCE} -q ${QUERY} -t ${THREADS} -n ${OUTPUT_NAME} -p ${PERCENT_ZEROES_FILTER} -o ${OUTPUT_DIR}