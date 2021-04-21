#!/bin/bash
#SBATCH --job-name=ecc_sv_finder_parallel_rev
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

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_sv_finder/guy11_v_all_AOC_repeatmasked



MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes/mapfile
COORDS_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/translocation_finder/guy11_v_all_repeatmasked/output_coords_rev
SIZES_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/translocation_finder/guy11_v_all_repeatmasked/output_sizes_rev
OUTPUT_DIR=/global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_sv_finder/guy11_v_all_AOC_repeatmasked/ecc_sv_finder_output_rev

if [ -d "${OUTPUT_DIR}" ]; then
    rm -r ${OUTPUT_DIR}
fi

mkdir ${OUTPUT_DIR}

if [ -f "jobqueue" ]; then
    rm "jobqueue"
fi

while read genome; do
    echo "python /global/home/users/pierrj/git/python/ecc_sv_finder_AOC.py ${COORDS_DIR}/${genome}.processed_output.coords ${SIZES_DIR}/${genome}.genomesize ${SIZES_DIR}/guy11_genome_baoetal2017.fasta.genomesize ${OUTPUT_DIR} ${genome} guy11" >> jobqueue
done < ${MAPFILE}

## change coords dir

parallel -j 24 < jobqueue

cd ${OUTPUT_DIR}

source activate mummer
module load gnuplot
module load imagemagick

REFERENCE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/guy11_genome_baoetal2017.fasta
GENOMES_PATH=/global/scratch/users/pierrj/eccDNA/pipeline_tests/sv_with_mummer/now_with_wheat_genomes_too/all_host_genomes
SUBDIR_OUTPUT=/global/scratch/users/pierrj/eccDNA/pipeline_tests/ecc_sv_finder/guy11_v_all_AOC_repeatmasked/ecc_sv_finder_output_rev/output_jpgs
PERCENT_ZEROES_FILTER=1

if [ -d "${SUBDIR_OUTPUT}" ]; then
    rm -r ${SUBDIR_OUTPUT}
fi

ls -1 | grep -v genome_mapfile | awk '{print substr($1,0,length($1)-8)}'> genome_mapfile

mkdir ${SUBDIR_OUTPUT}

while read reference;
do
## turn query name into path to genome
    reference_path=${GENOMES_PATH}/${reference}_genomic.fna
    cd ${reference}_v_guy11
    /global/home/users/pierrj/git/bash/mummerplotter_extract_coord.sh -r ${reference_path} -q ${REFERENCE} -e ${reference} -u guy11 -o ${SUBDIR_OUTPUT} -p ${PERCENT_ZEROES_FILTER}
    cd ..
done < genome_mapfile