#!/bin/bash
#SBATCH --job-name=gather_data_all_svs_wheat
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

WORKING_DIR=/global/scratch/users/pierrj/PAV_SV/SV/dels_wheat_blast_density_hist
PLOT_PROFILE_DIR=/global/scratch/users/pierrj/PAV_SV/SV/dels_wheat_blast
sv=DEL
DEL_BED=${PLOT_PROFILE_DIR}/all_${sv}.bed
SHUFFLED_BED=${PLOT_PROFILE_DIR}/shuffled_${sv}.bed
FASTA_FILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate/GCA_004785725.2_ASM478572v2_genomic.no_mito_no_mini.fna
GENOME_FILE=$FASTA_FILE
GENOME_CHROMSIZES=${FASTA_FILE}.chromsizes

cd $WORKING_DIR

# making prep files
for BED_FILE in $DEL_BED $SHUFFLED_BED; do
    if [[ "${BED_FILE}" == "${DEL_BED}" ]]; then
        label=del
    elif [[ "${BED_FILE}" == "${SHUFFLED_BED}" ]]; then
        label=shuffled
    fi
    awk '{print $3-$2}' ${BED_FILE} > ${label}.lengths
done

if [ -f "input_list" ]; then
    rm input_list
fi

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

# GC CONTENT
OUTPUT_NAME=gc_content
for BED_FILE in $DEL_BED $SHUFFLED_BED; do
    if [[ "${BED_FILE}" == "${DEL_BED}" ]]; then
        label=del
    elif [[ "${BED_FILE}" == "${SHUFFLED_BED}" ]]; then
        label=shuffled
    fi
    echo "bedtools nuc -fi ${GENOME_FILE} -bed $BED_FILE > ${label}.${OUTPUT_NAME}" >> jobqueue
done

# eccDNA seq
ECCDNA_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11_eccdnas/
BAMFILE_MAPFILE=${ECCDNA_DIR}/bamfile_mapfile
OUTPUT_NAME=eccdna
echo -e ${BAMFILE_MAPFILE}'\t'${OUTPUT_NAME} >> input_list

# in culture rnaseq
RNASEQ_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/
SRA_LIST=zhang_et_al_2019_complete_medium
BAMFILE_MAPFILE=${RNASEQ_DIR}/${SRA_LIST}_bamfile_mapfile
OUTPUT_NAME=inculture
echo -e ${BAMFILE_MAPFILE}'\t'${OUTPUT_NAME} >> input_list

# in planta rnaseq
RNASEQ_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/random_forest/guy11/
SRA_LIST=zhang_et_al_2019_in_planta
BAMFILE_MAPFILE=${RNASEQ_DIR}/${SRA_LIST}_bamfile_mapfile
OUTPUT_NAME=inplanta
echo -e ${BAMFILE_MAPFILE}'\t'${OUTPUT_NAME} >> input_list

## do all of the coverage calcs in parallel
while read BAMFILE_MAPFILE OUTPUT_NAME; do
    for BED_FILE in $DEL_BED $SHUFFLED_BED; do
        if [[ "${BED_FILE}" == "${DEL_BED}" ]]; then
            LABEL=del
        elif [[ "${BED_FILE}" == "${SHUFFLED_BED}" ]]; then
            LABEL=shuffled
        fi
        echo "/global/home/users/pierrj/git/bash/average_coverage_files.sh -b ${BED_FILE} -l ${LABEL} -m ${BAMFILE_MAPFILE} -g ${GENOME_CHROMSIZES} -o ${OUTPUT_NAME}" >> jobqueue
    done
done < input_list

parallel -j ${SLURM_NTASKS} < jobqueue