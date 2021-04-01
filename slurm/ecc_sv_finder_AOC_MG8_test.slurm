#!/bin/bash
#SBATCH --job-name=ecc_sv_finder_AOC_MG8_test
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_sv_finder/guy11_v_all_AOC_MG8_test/mummer_parse_output_all_false_positives

source activate mummer
module load gnuplot
module load imagemagick

REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.fasta
GENOME_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes
SUBDIR_OUTPUT=output_jpgs
PERCENT_ZEROES_FILTER=1

/global/home/users/pierrj/git/bash/auto_mummer_plot.sh -r ${REFERENCE} -g ${GENOME_DIR} -s ${SUBDIR_OUTPUT} -p ${PERCENT_ZEROES_FILTER}