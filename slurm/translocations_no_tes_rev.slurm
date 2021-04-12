#!/bin/bash
#SBATCH --job-name=translocations_no_tes_rev
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/translocation_finder/guy11_v_all_repeatmasked/

source ~/.bashrc
conda activate mummer
module load repeatmasker/4.1.1

MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/mapfile
REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/guy11_genome_baoetal2017.fasta
GENOMES_PATH=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/

OUTPUT_DIR_SIZES=output_sizes_rev
OUTPUT_DIR_BEDS=output_beds_rev
OUTPUT_DIR_COORDS=output_coords_rev
TMP_DIR=tmp_rev
RM_LIB=repbase_plus_guy11_modeler_rm_database.fasta

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
# # RepeatMasker -gff -nolow -no_is -norna -pa 24 -lib ${RM_LIB} -dir . ${REFERENCE}
# while read genome; do
#     samtools faidx ${GENOMES_PATH}/${genome}_genomic.fna
#     cut -f1,2 ${GENOMES_PATH}/${genome}_genomic.fna.fai > ${OUTPUT_DIR_SIZES}/$(basename ${genome}).genomesize
#     # RepeatMasker -gff -nolow -no_is -norna -pa 24 -lib ${RM_LIB} -dir . ${GENOMES_PATH}/${genome}_genomic.fna
#     /global/scratch/pierrj/mummer_4/bin/nucmer -t ${SLURM_NTASKS} --maxmatch -p ${TMP_DIR}/${genome} ${genome}_genomic.fna.masked $(basename ${REFERENCE}).masked
#     /global/scratch/pierrj/mummer_4/bin/show-coords ${TMP_DIR}/${genome}.delta > ${TMP_DIR}/${genome}.coords
#     tail -n+6 ${TMP_DIR}/${genome}.coords | awk -v OFS='\t' '{print $1, $2, $4, $5, $12, $13}' > ${OUTPUT_DIR_COORDS}/${genome}.processed_output.coords
# done < ${MAPFILE}

conda deactivate

COORDS_DIR=output_coords_rev
SIZES_DIR=output_sizes_rev
OUTPUT_DIR=translocation_finder_output_rev

if [ -d "${OUTPUT_DIR}" ]; then
    rm -r ${OUTPUT_DIR}
fi

mkdir ${OUTPUT_DIR}

if [ -f "jobqueue" ]; then
    rm "jobqueue"
fi

while read genome; do
    echo "python /global/home/users/pierrj/git/python/translocation_finder_AOC.py ${COORDS_DIR}/${genome}.processed_output.coords ${SIZES_DIR}/${genome}.genomesize ${SIZES_DIR}/guy11_genome_baoetal2017.fasta.genomesize ${OUTPUT_DIR} ${genome}  guy11" >> jobqueue
done < ${MAPFILE}

parallel -j 24 < jobqueue

cd ${OUTPUT_DIR}

conda activate mummer

module load gnuplot
module load imagemagick

REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.fasta
GENOME_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes
SUBDIR_OUTPUT=output_jpgs
PERCENT_ZEROES_FILTER=0.1

/global/home/users/pierrj/git/bash/auto_mummer_plot.sh -r ${REFERENCE} -g ${GENOME_DIR} -s ${SUBDIR_OUTPUT} -p ${PERCENT_ZEROES_FILTER}