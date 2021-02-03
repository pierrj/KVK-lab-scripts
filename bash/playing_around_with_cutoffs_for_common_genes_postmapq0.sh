#!/bin/bash
treatment=G3
MAPFILE=mapfile
GENE_BEDFILE=/global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.just_genes.post_repeatmasking.T0name.bed
OUTPUT_NAME=G3
PANNZER_INPUT=/global/scratch/users/pierrj/references/pannzer_GUY11_PacBio_merge_KROJ_post_repeatmasking.txt
SCORE_CUTOFF=0.6
TE_BEDFILE=/global/scratch/users/pierrj/references/GUY11_pacbio_ET_KROJ.features.justtranpsons.bed
SAMPLE_MAPFILE=G3.sample_mapfile
TREATMENT=G3

# while read SAMPLE; do cp /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/${SAMPLE}/ecccaller_output.${SAMPLE}.renamed.details.tsv . ; done < ${MAPFILE}

rm -r *

cp /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile .
cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/pva_comparison/G3.sample_mapfile .

while read sample; do
cp /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.confirmedsplitreads.bed .
done < mapfile

module remove r

module load r/3.4.2

if [ -f "${TREATMENT}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${TREATMENT}.mapfile_for_normalize_and_average_filecolumn
fi
while read SAMPLE; do
    bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b ${SAMPLE}.confirmedsplitreads.bed | awk -v OFS='\t' '{print $4, $5}' > ${SAMPLE}.splitreadspergene
    num_srs=$(wc -l ${SAMPLE}.confirmedsplitreads.bed | awk '{print $1/100000}')
    awk -v N=$num_srs '{print $1, $2/N}' ${SAMPLE}.splitreadspergene > ${SAMPLE}.normalized.splitreadspergene ## NORMALIZE TO DEAL WITH FAVORING OF SMALLER GENES TEST THIS LATER
    echo ${SAMPLE}.normalized.splitreadspergene >> ${TREATMENT}.mapfile_for_normalize_and_average_filecolumn
done < ${MAPFILE}

if [ -f "${TREATMENT}.normalize_table_column" ]; then
    rm ${TREATMENT}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${TREATMENT}.normalize_table_column ; done
paste ${TREATMENT}.mapfile_for_normalize_and_average_filecolumn ${TREATMENT}.normalize_table_column ${SAMPLE_MAPFILE} > ${TREATMENT}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${TREATMENT}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${TREATMENT}.normalized_binned ${TREATMENT}.normalized.splitreadspergene

sort -k2,2nr ${TREATMENT}.normalized.splitreadspergene | head -1000 | awk '{print $1}' > ${TREATMENT}.common.genes

awk '$3==0' ${TREATMENT}.normalized.splitreadspergene > ${TREATMENT}.neverfound.genes

if [ -d "raw_files" ]; then
    rm -r raw_files
fi

mkdir raw_files

mv * raw_files/ 2>/dev/null

if [ -d "common_distance" ]; then
    rm -r common_distance
fi

mkdir common_distance

cd common_distance

Rscript --vanilla /global/home/users/pierrj/git/R/nearest_neighbor_calc_forpipeline.R ${GENE_BEDFILE} ${GENE_BEDFILE} ../raw_files/${OUTPUT_NAME}.common.genes ${OUTPUT_NAME}.common_genes.genedistance.distanceplot.pdf ${OUTPUT_NAME}.common_genes.genedistance.pvalues.txt

Rscript --vanilla /global/home/users/pierrj/git/R/nearest_neighbor_calc_forpipeline.R ${GENE_BEDFILE} ${TE_BEDFILE} ../raw_files/${OUTPUT_NAME}.common.genes ${OUTPUT_NAME}.common_genes.tedistance.distanceplot.pdf ${OUTPUT_NAME}.common_genes.tedistance.pvalues.txt

cd ..

if [ -d "never_distance" ]; then
    rm -r never_distance
fi

mkdir never_distance

cd never_distance

Rscript --vanilla /global/home/users/pierrj/git/R/nearest_neighbor_calc_forpipeline.R ${GENE_BEDFILE} ${GENE_BEDFILE} ../raw_files/${OUTPUT_NAME}.neverfound.genes ${OUTPUT_NAME}.neverfound_genes.genedistance.distanceplot.pdf ${OUTPUT_NAME}.neverfound_genes.genedistance.pvalues.txt

Rscript --vanilla /global/home/users/pierrj/git/R/nearest_neighbor_calc_forpipeline.R ${GENE_BEDFILE} ${TE_BEDFILE} ../raw_files/${OUTPUT_NAME}.neverfound.genes ${OUTPUT_NAME}.neverfound_genes.tedistance.distanceplot.pdf ${OUTPUT_NAME}.neverfound_genes.tedistance.pvalues.txt

cd ..

if [ -d "common_goenrichment" ]; then
    rm -r common_goenrichment
fi

mkdir common_goenrichment

cd common_goenrichment

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.common_genes.goenrichment.${SCORE_CUTOFF}_MF -c ${SCORE_CUTOFF} -s ../raw_files/${OUTPUT_NAME}.common.genes -g MF

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.common_genes.goenrichment.${SCORE_CUTOFF}_BP -c ${SCORE_CUTOFF} -s ../raw_files/${OUTPUT_NAME}.common.genes -g BP

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.common_genes.goenrichment.${SCORE_CUTOFF}_CC -c ${SCORE_CUTOFF} -s ../raw_files/${OUTPUT_NAME}.common.genes -g CC

cd ..

if [ -d "never_goenrichment" ]; then
    rm -r never_goenrichment
fi

mkdir never_goenrichment

cd never_goenrichment

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.neverfound_genes.goenrichment.${SCORE_CUTOFF}_MF -c ${SCORE_CUTOFF} -s ../raw_files/${OUTPUT_NAME}.neverfound.genes -g MF

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.neverfound_genes.goenrichment.${SCORE_CUTOFF}_BP -c ${SCORE_CUTOFF} -s ../raw_files/${OUTPUT_NAME}.neverfound.genes -g BP

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.neverfound_genes.goenrichment.${SCORE_CUTOFF}_CC -c ${SCORE_CUTOFF} -s ../raw_files/${OUTPUT_NAME}.neverfound.genes -g CC

cd ..