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

source /global/scratch/users/pierrj/spaced_anchor/mugsy/mugsyenv.sh

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

if [ -d "${TMP_DIR}" ]; then
    rm -r ${TMP_DIR}
fi
mkdir ${TMP_DIR}

if [ ! -f "${REFERENCE}.fai" ]; then
    samtools faidx ${REFERENCE}
fi
cut -f1,2 ${REFERENCE}.fai > ${TMP_DIR}/${ref}.genomesize
if [ ! -f "${QUERY}.fai" ]; then
    samtools faidx ${QUERY}
fi
cut -f1,2 ${QUERY}.fai > ${TMP_DIR}/${quer}.genomesize
mugsy --directory . --prefix ${OUTPUT_NAME} ${REFERENCE} ${QUERY}
/global/scratch/pierrj/mummer_4/bin/show-coords ${OUTPUT_NAME}.guy11_genome_baoetal2017.delta > ${TMP_DIR}/${OUTPUT_NAME}.coords
tail -n+6 ${OUTPUT_NAME}.coords | awk -v OFS='\t' '{print $1, $2, $4, $5, $12, $13}' > ${TMP_DIR}/${OUTPUT_NAME}.processed_output.coords

python /global/home/users/pierrj/git/python/ecc_sv_finder_AOC.py ${TMP_DIR}/${OUTPUT_NAME}.processed_output.coords ${TMP_DIR}/${ref}.genomesize ${TMP_DIR}/${quer}.genomesize ${TMP_DIR} ${ref} ${quer}