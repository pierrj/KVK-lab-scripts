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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/test_ecc_sv_finder_script/magor_neurospora_test

REFERENCE=guy11_genome_baoetal2017.fasta
QUERY=Neucr2_AssemblyScaffolds.fasta
THREADS=${SLURM_NTASKS}
OUTPUT_NAME=guy11_v_neuro
PERCENT_ZEROES_FILTER=1.0
OUTPUT_DIR=guy11_v_neuro_jpg
TMP_DIR=guy11_v_neuro_tmp

source ~/.bashrc
conda activate mummer

module load gnuplot
module load imagemagick

/global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r ${REFERENCE} -q ${QUERY} -t ${THREADS} -n ${OUTPUT_NAME} -p ${PERCENT_ZEROES_FILTER} -o ${OUTPUT_DIR} -m ${TMP_DIR}


# REFERENCE=M22_yeast.fasta
# QUERY=EC1118_yeast.fasta
# THREADS=${SLURM_NTASKS}
# OUTPUT_NAME=22_v_1118
# PERCENT_ZEROES_FILTER=1.0
# OUTPUT_DIR=22_v_1118_jpg
# TMP_DIR=22_v_1118_tmp

# /global/home/users/pierrj/git/bash/ecc_sv_finder.sh -r ${REFERENCE} -q ${QUERY} -t ${THREADS} -n ${OUTPUT_NAME} -p ${PERCENT_ZEROES_FILTER} -o ${OUTPUT_DIR} -m ${TMP_DIR}