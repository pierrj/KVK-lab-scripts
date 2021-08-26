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

# cp /global/scratch/pierrj/eccDNA/rerunning_stuff_final/comparative/count_comparison/ecc_count_metadata_unix.txt .

# cp /global/scratch/pierrj/eccDNA/rerunning_stuff_final/comparative/count_comparison/mapfile_ecc_count .

# uniq ecc_count_metadata_unix.txt > ecc_count_metadata_unix.txt.uniq

# paste <(awk -v OFS='=' '{print $1, $2, $3, $4, $5}' ecc_count_metadata_unix.txt.uniq) <(awk -v OFS='\t' '{print $3, $6, $1, $2, $4, $5}' mapfile_ecc_count) > mapfile_gene_portions

module load emboss

source activate agat

export PERL5LIB=~/.conda/envs/agat/lib

while read -r ORGANISM_NAME GENE_GFF SAMPLE_MAPFILE ECC_DIR GENOME_FILE COPIA_FILE GYPSY_FILE CONTIGNAMES_FILE
do
    if [ -d "${ORGANISM_NAME}" ]; then
        rm -r ${ORGANISM_NAME}
    fi
    mkdir ${ORGANISM_NAME}
    cd ${ORGANISM_NAME}
        /global/home/users/pierrj/git/bash/gene_portions_enrichment.sh -o ${ORGANISM_NAME} \
            -g ${GENE_GFF} -s ${SAMPLE_MAPFILE} -d ${ECC_DIR} \
            -f ${GENOME_FILE} -c ${COPIA_FILE} -y ${GYPSY_FILE} -n ${CONTIGNAMES_FILE}
    cd ..
done < mapfile_gene_portions

if [ -f "gene_portions_full_table" ]; then
    rm gene_portions_full_table
fi

while read sample; do
cd ${sample}
    if [ -f "temp.sample_name_column" ]; then
        rm temp.sample_name_column
    fi
    for i in {0..8}; do
        echo ${sample} >> temp.sample_name_column
    done
    awk '{print $1}' observed_averages_micro > temp.gene_portion_name_column
    awk '{print $2}' observed_averages_micro > temp.observed_averages_micro_column
    awk '{print $2}' expected_averages_micro > temp.expected_averages_micro_column
    awk '{print $2}' observed_averages_large > temp.observed_averages_large_column
    awk '{print $2}' expected_averages_large > temp.expected_averages_large_column
    paste temp.sample_name_column temp.gene_portion_name_column temp.observed_averages_micro_column temp.expected_averages_micro_column temp.observed_averages_large_column temp.expected_averages_large_column >> ../gene_portions_full_table
    rm temp.*
cd ..
done < mapfile_studies