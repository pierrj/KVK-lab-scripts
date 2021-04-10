#!/bin/bash
#SBATCH --job-name=translocation_finder_parallel
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/translocation_finder/guy11_v_all

MAPFILE=/global/scratch/pierrj/eccDNA/pipeline_tests/translocation_finder/guy11_v_all/new_new_mapfile
COORDS_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/translocation_finder/guy11_v_all/output_coords
SIZES_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/translocation_finder/guy11_v_all/output_sizes
OUTPUT_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/translocation_finder/guy11_v_all/translocation_finder_output

if [ -d "${OUTPUT_DIR}" ]; then
    rm -r ${OUTPUT_DIR}
fi

mkdir ${OUTPUT_DIR}

if [ -f "jobqueue" ]; then
    rm "jobqueue"
fi

while read genome; do
    echo "python /global/home/users/pierrj/git/python/translocation_finder.py ${COORDS_DIR}/${genome}.processed_output.coords ${SIZES_DIR}/guy11_genome_baoetal2017.fasta.genomesize ${SIZES_DIR}/${genome}.genomesize ${OUTPUT_DIR} guy11 ${genome}" >> jobqueue
done < ${MAPFILE}

parallel -j 24 < jobqueue

cd ${OUTPUT_DIR}

source activate mummer
module load gnuplot
module load imagemagick

REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.fasta
GENOME_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes
SUBDIR_OUTPUT=output_jpgs
PERCENT_ZEROES_FILTER=0.1

/global/home/users/pierrj/git/bash/auto_mummer_plot.sh -r ${REFERENCE} -g ${GENOME_DIR} -s ${SUBDIR_OUTPUT} -p ${PERCENT_ZEROES_FILTER}