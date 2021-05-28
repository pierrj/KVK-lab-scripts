#!/bin/bash
#SBATCH --job-name=gene_analysis_final
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis

source activate pfamscan

module remove r

module load r/3.4.2

PANNZER_INPUT=/global/scratch/users/pierrj/references/pannzer_output_fungap_guy11_12_28_20.txt
OUTPUT_NAME=G3
SCORE_CUTOFF=0.6
MAPFILE=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile
GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.renamed.bed
TE_BEDFILE=
REPEAT_BEDFILE=
PFAM_DIR=/global/scratch/users/pierrj/pfam/
CDS_FASTA=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20_prot.faa
THREADS=${SLURM_NTASKS}
SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_comparison/G3.sample_mapfile

if [ -f "${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn
fi
while read SAMPLE; do
    bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b ${SAMPLE}.confirmedsplitreads.bed | awk -v OFS='\t' '{print $4, $5}' > ${SAMPLE}.splitreadspergene
    num_srs=$(wc -l ${SAMPLE}.confirmedsplitreads.bed | awk '{print $1/100000}')
    awk -v N=$num_srs '{print $1, $2/N}' ${SAMPLE}.splitreadspergene > ${SAMPLE}.normalized.splitreadspergene ## NORMALIZE TO DEAL WITH FAVORING OF SMALLER GENES TEST THIS LATER
    echo ${SAMPLE}.normalized.splitreadspergene >> ${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn
done < ${MAPFILE}

if [ -f "${OUTPUT_NAME}.normalize_table_column" ]; then
    rm ${OUTPUT_NAME}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${OUTPUT_NAME}.normalize_table_column ; done
paste ${OUTPUT_NAME}.mapfile_for_normalize_and_average_filecolumn ${OUTPUT_NAME}.normalize_table_column ${SAMPLE_MAPFILE} > ${OUTPUT_NAME}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${OUTPUT_NAME}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${OUTPUT_NAME}.normalized_binned ${OUTPUT_NAME}.normalized.splitreadspergene

top_ten_percent=$(awk '$3!=0' ${OUTPUT_NAME}.normalized.splitreadspergene | wc -l | awk '{print int($1/10)}')
sort -k2,2nr ${OUTPUT_NAME}.normalized.splitreadspergene | head -${top_ten_percent} | awk '{print $1}' > common/${OUTPUT_NAME}.common.genes

awk '{ if ($3==0) {print $1}}' ${OUTPUT_NAME}.normalized.splitreadspergene > neverfound/${OUTPUT_NAME}.neverfound.genes

awk '{print $1}' ${OUTPUT_NAME}.normalized.splitreadspergene > common/${OUTPUT_NAME}.allgenenames

awk '{print $1}' ${OUTPUT_NAME}.normalized.splitreadspergene > neverfound/${OUTPUT_NAME}.allgenenames

cp ${GENE_BEDFILE} common/
cp ${GENE_BEDFILE} neverfound/

cp ${TE_BEDFILE} common/
cp ${TE_BEDFILE} neverfound/

cp ${REPEAT_BEDFILE} common/
cp ${REPEAT_BEDFILE} neverfound/

cp ${PANNZER_INPUT} common/
cp ${PANNZER_INPUT} neverfound/


python /global/home/users/pierrj/git/python/parse_pannzer_output.py ${PANNZER_INPUT} ${SCORE_CUTOFF} ${OUTPUT_NAME}.filtered.GO.txt

awk -v OFS='\t' '{print $0}' ${OUTPUT_NAME}.filtered.GO.txt > ${OUTPUT_NAME}.filtered.GO.awk.txt

cp ${OUTPUT_NAME}.filtered.GO.awk.txt common/
cp ${OUTPUT_NAME}.filtered.GO.awk.txt neverfound/


cd common

/global/home/users/pierrj/git/bash/pfam_scan_wordcloud_wrapper.sh -g ${OUTPUT_NAME}.common.genes -a ${OUTPUT_NAME}.allgenenames -n ${OUTPUT_NAME}.common_genes -p ${PFAM_DIR} -c ${CDS_FASTA} -t ${THREADS}

cd ..


cd neverfound

/global/home/users/pierrj/git/bash/pfam_scan_wordcloud_wrapper.sh -g ${OUTPUT_NAME}.neverfound.genes -a ${OUTPUT_NAME}.allgenenames -n ${OUTPUT_NAME}.never_genes -p ${PFAM_DIR} -c ${CDS_FASTA} -t ${THREADS}

cd ..