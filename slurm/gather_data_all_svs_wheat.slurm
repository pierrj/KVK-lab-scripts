#!/bin/bash
#SBATCH --job-name=gather_data_all_svs_wheat
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

WORKING_DIR=/global/scratch/users/pierrj/PAV_SV/SV/dels_wheat_blast_density_hist
PLOT_PROFILE_DIR=/global/scratch/users/pierrj/PAV_SV/SV/dels_wheat_blast
sv=DEL
DEL_BED=${PLOT_PROFILE_DIR}/all_${sv}.bed
FASTA_FILE=/global/scratch/users/pierrj/fungap_runs/wheat_blast/genomes_to_annotate/GCA_004785725.2_ASM478572v2_genomic.no_mito_no_mini.fna
GENOME_FILE=$FASTA_FILE
GENOME_CHROMSIZES=${FASTA_FILE}.chromsizes

cd $WORKING_DIR

# making prep files

## only keep the ~90% shortest deletions and sort

DEL_BED_ORIGINAL=$DEL_BED

awk '{ if ($3-$2 < 10000) {print $0}}' $DEL_BED | sort -k1,1 -k2,2n > dels.subset.bed
DEL_BED=dels.subset.bed

## exclude ALL original deletions but shuffle only the small deletions
bedtools shuffle -i $DEL_BED -g $GENOME_CHROMSIZES -excl $DEL_BED_ORIGINAL | sort -k1,1 -k2,2n > shuffled.bed
SHUFFLED_BED=shuffled.bed

for BED_FILE in $DEL_BED $SHUFFLED_BED; do
    if [[ "${BED_FILE}" == "${DEL_BED}" ]]; then
        LABEL=del
    elif [[ "${BED_FILE}" == "${SHUFFLED_BED}" ]]; then
        LABEL=shuffled
    fi
    awk '{print $3-$2}' ${BED_FILE} > ${LABEL}.lengths
done

if [ -f "input_list" ]; then
    rm input_list
fi

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

# in culture rnaseq
RNASEQ_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/random_forest/b71
SRA_LIST=b71_assembly_complete_medium
BAMFILE_MAPFILE=${RNASEQ_DIR}/${SRA_LIST}_bamfile_mapfile
OUTPUT_NAME=inculture
echo -e ${BAMFILE_MAPFILE}'\t'${OUTPUT_NAME} >> input_list

# in planta rnaseq
RNASEQ_DIR=/global/scratch/users/pierrj/PAV_SV/PAV/wheat_blast_all/random_forest/b71
SRA_LIST=b71_assembly_complete_medium_in_planta
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

## some cleaning before starting the job just in case
rm *coverage*

parallel -j ${SLURM_NTASKS} < jobqueue

# GC CONTENT using bedtools nuc
OUTPUT_NAME=gc_content
for BED_FILE in $DEL_BED $SHUFFLED_BED; do
    if [[ "${BED_FILE}" == "${DEL_BED}" ]]; then
        LABEL=del
    elif [[ "${BED_FILE}" == "${SHUFFLED_BED}" ]]; then
        LABEL=shuffled
    fi
    ### THIS IS ACTUALLY AT CONTENT NOT GC CONTENT ###
    bedtools nuc -fi ${GENOME_FILE} -bed $BED_FILE | awk -v OFS='\t' '{ if (NR > 1) {print $4}}' > ${LABEL}.${OUTPUT_NAME}.txt
done