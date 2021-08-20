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

rm org*

rm mapfile*

# paste <(awk -v OFS='\t' '{print $3, $6, $1, $2, $4, $5}' mapfile_ecc_count) org_names > mapfile_gene_portions

# module load genometools
# module load emboss

# source activate agat

# export PERL5LIB=~/.conda/envs/agat/lib



# while read -r dir GENE_GFF SAMPLE_MAPFILE ECC_DIR GENOME_FILE COPIA_FILE GYPSY_FILE ORGANISM_NAME
# do
#     if [ -d "${ORGANISM_NAME}" ]; then
#         rm -r ${ORGANISM_NAME}
#     fi
#     mkdir ${ORGANISM_NAME}
#     cd ${ORGANISM_NAME}
#     echo ${ORGANISM_NAME}
#     echo ${GENE_GFF}
#     echo ${SAMPLE_MAPFILE}
#     echo ${ECC_DIR}
#     echo ${GENOME_FILE}
#     echo ${COPIA_FILE}
#     echo ${GYPSY_FILE}
#         # /global/home/users/pierrj/git/bash/gene_portions_enrichment.sh -o ${ORGANISM_NAME} \
#         #     -g ${GENE_GFF} -s ${SAMPLE_MAPFILE} -d ${ECC_DIR} \
#         #     -f ${GENOME_FILE} -c ${COPIA_FILE} -y ${GYPSY_FILE}
#     cd ..
# done < mapfile_gene_portions