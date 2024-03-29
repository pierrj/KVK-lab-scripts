#!/bin/bash
#SBATCH --job-name=mummer_guy11_every_genome
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too

source activate mummer

MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/mapfile
REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/guy11_genome_baoetal2017.fasta
OUTPUT_DIR_SIZES=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/output_sizes_rev
OUTPUT_DIR_BEDS=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/output_beds_rev
OUTPUT_DIR_COORDS=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/output_coords_rev
TMP_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/tmp_rev
GENOMES_PATH=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/

if [ -d "${OUTPUT_DIR_SIZES}" ]; then
    rm -r ${OUTPUT_DIR_SIZES}
fi

if [ -d "${OUTPUT_DIR_BEDS}" ]; then
    rm -r ${OUTPUT_DIR_BEDS}
fi

if [ -d "${OUTPUT_DIR_COORDS}" ]; then
    rm -r ${OUTPUT_DIR_COORDS}
fi

if [ -d "${TMP_DIR}" ]; then
    rm -r ${TMP_DIR}
fi

mkdir ${OUTPUT_DIR_SIZES}
mkdir ${OUTPUT_DIR_BEDS}
mkdir ${TMP_DIR}
mkdir ${OUTPUT_DIR_COORDS}

samtools faidx ${REFERENCE}
cut -f1,2 ${REFERENCE}.fai > ${OUTPUT_DIR_SIZES}/$(basename ${REFERENCE}).genomesize
while read genome; do
    samtools faidx ${GENOMES_PATH}/${genome}_genomic.fna
    cut -f1,2 ${GENOMES_PATH}/${genome}_genomic.fna.fai > ${OUTPUT_DIR_SIZES}/$(basename ${genome}).genomesize
    /global/scratch/pierrj/mummer_4/bin/nucmer -t ${SLURM_NTASKS} --maxmatch -p ${TMP_DIR}/${genome} ${GENOMES_PATH}/${genome}_genomic.fna ${GENOMES_PATH}/guy11_genome_baoetal2017.fasta
    /global/scratch/pierrj/mummer_4/bin/show-coords ${TMP_DIR}/${genome}.delta > ${TMP_DIR}/${genome}.coords
    tail -n+6 ${TMP_DIR}/${genome}.coords | awk -v OFS='\t' '{print $1, $2, $4, $5, $12, $13}' > ${OUTPUT_DIR_COORDS}/${genome}.processed_output.coords
done < ${MAPFILE}

# OUTPUT_DIR_SIZES=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/output_sizes
# OUTPUT_DIR_BEDS=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/output_beds
# OUTPUT_DIR_COORDS=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/output_coords
# TMP_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/tmp

# if [ -d "${OUTPUT_DIR_SIZES}" ]; then
#     rm -r ${OUTPUT_DIR_SIZES}
# fi

# if [ -d "${OUTPUT_DIR_BEDS}" ]; then
#     rm -r ${OUTPUT_DIR_BEDS}
# fi

# if [ -d "${OUTPUT_DIR_COORDS}" ]; then
#     rm -r ${OUTPUT_DIR_COORDS}
# fi

# if [ -d "${TMP_DIR}" ]; then
#     rm -r ${TMP_DIR}
# fi

# mkdir ${OUTPUT_DIR_SIZES}
# mkdir ${OUTPUT_DIR_BEDS}
# mkdir ${TMP_DIR}
# mkdir ${OUTPUT_DIR_COORDS}

# samtools faidx ${REFERENCE}
# cut -f1,2 ${REFERENCE}.fai > ${OUTPUT_DIR_SIZES}/$(basename ${REFERENCE}).genomesize
# while read genome; do
#     samtools faidx ${GENOMES_PATH}/${genome}_genomic.fna
#     cut -f1,2 ${GENOMES_PATH}/${genome}_genomic.fna.fai > ${OUTPUT_DIR_SIZES}/$(basename ${genome}).genomesize
#     /global/scratch/pierrj/mummer_4/bin/nucmer -t ${SLURM_NTASKS} --maxmatch -p ${TMP_DIR}/${genome} ${GENOMES_PATH}/guy11_genome_baoetal2017.fasta ${GENOMES_PATH}/${genome}_genomic.fna
#     /global/scratch/pierrj/mummer_4/bin/show-coords ${TMP_DIR}/${genome}.delta > ${TMP_DIR}/${genome}.coords
#     tail -n+6 ${TMP_DIR}/${genome}.coords | awk -v OFS='\t' '{print $1, $2, $4, $5, $12, $13}' > ${OUTPUT_DIR_COORDS}/${genome}.processed_output.coords
# done < ${MAPFILE}