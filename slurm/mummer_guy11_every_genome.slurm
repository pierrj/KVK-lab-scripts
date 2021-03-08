#!/bin/bash
#SBATCH --job-name=mummer_guy11_every_genome
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too

cd all_host_genomes
gunzip *
cd ..

# source activate mummer

# MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/moryzae_genomes/mapfile
# REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/moryzae_genomes/guy11_genome_baoetal2017.fasta
# GENOMES_PATH=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/moryzae_genomes/
# OUTPUT_DIR_SIZES=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/output_sizes
# OUTPUT_DIR_BEDS=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/output_beds
# TMP_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/tmp

# if [ -d "${OUTPUT_DIR}" ]; then
#     rm -r ${OUTPUT_DIR}
# fi

# if [ -d "${TMP_DIR}" ]; then
#     rm -r ${TMP_DIR}
# fi

# mkdir ${OUTPUT_DIR}
# mkdir ${TMP_DIR}

# samtools faidx ${REFERENCE}
# cut -f1,2 ${REFERENCE}.fai > ${OUTPUT_DIR}/$(basename ${REFERENCE}).genomesize
# while read genome; do
#     samtools faidx ${GENOMES_PATH}/${genome}_genomic.fna
#     cut -f1,2 ${GENOMES_PATH}/${genome}_genomic.fna.fai > ${OUTPUT_DIR_SIZES}/$(basename ${genome}).genomesize
#     nucmer -p ${TMP_DIR}/${genome} ${GENOMES_PATH}/guy11_genome_baoetal2017.fasta ${GENOMES_PATH}/${genome}_genomic.fna
#     show-coords ${TMP_DIR}/${genome}.delta > ${TMP_DIR}/${genome}.coords
#     tail -n +6 ${TMP_DIR}/${genome}.coords | awk -v OFS='\t' '{print $12, $1, $2}' | sort -k1,1 -k2,2n >  ${TMP_DIR}/${genome}.bed
#     bedtools genomecov -bga -i ${TMP_DIR}/${genome}.bed -g ${OUTPUT_DIR}/$(basename ${REFERENCE}).genomesize | awk '$2!=0' | awk -v OFS='\t' '{if ($4==1) { print $1, $2, $3, $3-$2}}' > ${TMP_DIR}/single_alignment_regions.${genome}.bed
#     tail -n +6 ${TMP_DIR}/${genome}.coords | awk -v OFS='\t' '{print $12, $1, $2, $13, $4, $5}' | sort -k1,1 -k2,2n > ${TMP_DIR}/${genome}.alignments.bed
#     bedtools intersect -f 1 -wa -wb -a ${TMP_DIR}/single_alignment_regions.${genome}.bed -b ${TMP_DIR}/${genome}.alignments.bed > ${TMP_DIR}/${genome}.single_coverage.bed
#     bedtools genomecov -bga -i ${TMP_DIR}/${genome}.bed -g ${OUTPUT_DIR_SIZES}/$(basename ${REFERENCE}).genomesize | awk '$2!=0' | awk -v OFS='\t' '{if ($4==0) { print $1, $2, $3,  $3-$2, "no_match"}}' > ${TMP_DIR}/${genome}.no_match.bed
#     bedtools genomecov -bga -i ${TMP_DIR}/${genome}.bed -g ${OUTPUT_DIR_SIZES}/$(basename ${REFERENCE}).genomesize | awk '$2!=0' | awk -v OFS='\t' '{if ($4>1) { print $1, $2, $3,  $3-$2, "repetitive"}}' > ${TMP_DIR}/${genome}.repetitive.bed
#     cat ${TMP_DIR}/${genome}.repetitive.bed ${TMP_DIR}/${genome}.no_match.bed ${TMP_DIR}/${genome}.single_coverage.bed | sort -k1,1 -k2,2n > ${OUTPUT_DIR_BEDS}/${genome}.processed_output.bed
# done < ${MAPFILE}