#!/bin/bash
#SBATCH --job-name=gene_portions_all
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/script_test

module load genometools
module load emboss

source activate agat

export PERL5LIB=~/.conda/envs/agat/lib

ORGANISM_NAME=human
GENE_GFF=/global/scratch/users/pierrj/references/GRCh37.p13.gff
SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/2018_moller/full_run/mapfile_muscle
ECC_DIR=/global/scratch/users/pierrj/eccDNA/2018_moller/full_run/
GENOME_FILE=/global/scratch/users/pierrj/references/GRCh37.p13.fasta
COPIA_FILE=/global/scratch/pierrj/references/te_annotations/human/blank
GYPSY_FILE=/global/scratch/pierrj/references/te_annotations/human/blank

if [ -d "${ORGANISM_NAME}" ]; then
    rm -r ${ORGANISM_NAME}
fi

mkdir ${ORGANISM_NAME}
cd ${ORGANISM_NAME}
    /global/home/users/pierrj/git/bash/gene_portions_enrichment.sh -o ${ORGANISM_NAME} \
        -g ${GENE_GFF} -s ${SAMPLE_MAPFILE} -d ${ECC_DIR} \
        -f ${GENOME_FILE} -c ${COPIA_FILE} -y ${GYPSY_FILE}
cd ..