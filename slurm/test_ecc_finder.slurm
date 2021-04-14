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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/test_ecc_sv_finder_script/yeast_test

REFERENCE=EC1118_yeast.fasta
QUERY=M22_yeast.fasta
THREADS=${SLURM_NTASKS}
OUTPUT_NAME=1118_v_22
PERCENT_ZEROES_FILTER=1.0
OUTPUT_DIR=118_v_22_jpg
TMP_DIR=1118_v_22_tmp

source ~/.bashrc
conda activate mummer

module load gnuplot
module load imagemagick

/global/home/users/pierrj/git/bash/ecc_sv_finder_maxmatch.sh -r ${REFERENCE} -q ${QUERY} -t ${THREADS} -n ${OUTPUT_NAME} -p ${PERCENT_ZEROES_FILTER} -o ${OUTPUT_DIR} -m ${TMP_DIR}


REFERENCE=M22_yeast.fasta
QUERY=EC1118_yeast.fasta
THREADS=${SLURM_NTASKS}
OUTPUT_NAME=22_v_1118
PERCENT_ZEROES_FILTER=1.0
OUTPUT_DIR=22_v_1118_jpg
TMP_DIR=22_v_1118_tmp

/global/home/users/pierrj/git/bash/ecc_sv_finder_maxmatch.sh -r ${REFERENCE} -q ${QUERY} -t ${THREADS} -n ${OUTPUT_NAME} -p ${PERCENT_ZEROES_FILTER} -o ${OUTPUT_DIR} -m ${TMP_DIR}