#!/bin/bash
while getopts p:n:c:m:b:e:f:a:t:i: option
do
case "${option}"
in
p) PANNZER_INPUT=${OPTARG};;
n) OUTPUT_NAME=${OPTARG};;
c) SCORE_CUTOFF=${OPTARG};;
m) MAPFILE=${OPTARG};;
b) GENE_BEDFILE=${OPTARG};;
e) TE_BEDFILE=${OPTARG};;
f) PFAM_DIR=${OPTARG};;
a) CDS_FASTA=${OPTARG};;
t) THREADS=${OPTARG};;
i) SAMPLE_MAPFILE=${OPTARG};;
esac
done

echo "THE START OF THIS SCRIPT IS OUTDATED PLEASE DONT USE"
echo "THE START OF THIS SCRIPT IS OUTDATED PLEASE DONT USE"

echo "THE START OF THIS SCRIPT IS OUTDATED PLEASE DONT USE"

echo "THE START OF THIS SCRIPT IS OUTDATED PLEASE DONT USE"



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

## maybe something besides 1000 here?
# maybe top 10% of found genes like this
top_ten_percent=$(awk '$3!=0' ${OUTPUT_NAME}.normalized.splitreadspergene | wc -l | awk '{print int($1/10)}')
sort -k2,2nr ${OUTPUT_NAME}.normalized.splitreadspergene | head -${top_ten_percent} | awk '{print $1}' > ${OUTPUT_NAME}.common.genes

awk '{ if ($3==0) {print $1}}' ${OUTPUT_NAME}.normalized.splitreadspergene > ${OUTPUT_NAME}.neverfound.genes

awk '{print $1}' ${OUTPUT_NAME}.normalized.splitreadspergene > ${OUTPUT_NAME}.allgenenames

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

if [ -d "common_wordcloud" ]; then
    rm -r common_wordcloud
fi

mkdir common_wordcloud

cd common_wordcloud

/global/home/users/pierrj/git/bash/pfam_scan_wordcloud_wrapper.sh -g ../raw_files/${OUTPUT_NAME}.common.genes -a ../raw_files/${OUTPUT_NAME}.allgenenames -n ${OUTPUT_NAME}.common_genes -p ${PFAM_DIR} -c ${CDS_FASTA} -t ${THREADS}

cd ..

if [ -d "never_wordcloud" ]; then
    rm -r never_wordcloud
fi

mkdir never_wordcloud

cd never_wordcloud

/global/home/users/pierrj/git/bash/pfam_scan_wordcloud_wrapper.sh -g ../raw_files/${OUTPUT_NAME}.neverfound.genes -a ../raw_files/${OUTPUT_NAME}.allgenenames -n ${OUTPUT_NAME}.never_genes -p ${PFAM_DIR} -c ${CDS_FASTA} -t ${THREADS}

cd ..