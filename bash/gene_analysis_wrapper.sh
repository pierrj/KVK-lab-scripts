#!/bin/bash
while getopts p:n:c:m:b:e:f:a:e: option
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
esac
done

while read SAMPLE
do
    awk '$6=="hconf" && $3-$2>1000' ecccaller_output.${SAMPLE}.renamed.details.tsv > ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv
    line_count=$(wc -l ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv | awk '{print int($1/100)}') 
    sr_count_cutoff=$(sort -k4,4nr ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv | head -${line_count} | tail -1 | awk '{print $4}')
    awk -v CUTOFF=${sr_count_cutoff} '$4>=CUTOFF' ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv > ecccaller_output.${SAMPLE}.common.renamed.details.tsv
done < ${MAPFILE}

## TO GENERATE BIOREP MAPFILE IF STANDARD FILES ARE SET

awk '{print substr($1, 0,4)}' ${MAPFILE} | sort | uniq > tmp_biorepmapfile

while read bio_rep; 
do
    cat $(find . -maxdepth 1 -name "ecccaller_output.${bio_rep}*.common.renamed.details.tsv" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") > ecccaller_output.${bio_rep}.common.renamed.details.tsv
done < tmp_biorepmapfile

while read bio_rep; 
do
    bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b ecccaller_output.${bio_rep}.common.renamed.details.tsv | awk '{if ($5!=0) {print $4}}' > ${bio_rep}.common.genes
done < tmp_biorepmapfile

if [ -f "${OUTPUT_NAME}.common.genes" ]; then
    rm ${OUTPUT_NAME}.common.genes
fi

cat $(find . -maxdepth 1 -name "*.common.genes" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | sort | uniq -c | awk '$1==3 {print $2}' > ${OUTPUT_NAME}.common.genes

while read bio_rep; 
do
    cat $(find . -maxdepth 1 -name "ecccaller_output.${bio_rep}*.hconf.renamed.details.tsv" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") > ecccaller_output.${bio_rep}.hconf.renamed.details.tsv
done < tmp_biorepmapfile

while read bio_rep; 
do
    bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b ecccaller_output.${bio_rep}.hconf.renamed.details.tsv | awk '{if ($5!=0) {print $4}}' > ${bio_rep}.hconf.genes
done < tmp_biorepmapfile


if [ -f "${OUTPUT_NAME}.hconf.genes" ]; then
    rm ${OUTPUT_NAME}.hconf.genes
fi

cat $(find . -maxdepth 1 -name "*.hconf.genes" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | sort | uniq > ${OUTPUT_NAME}.allhconf.genes

awk '{print $4}' ${GENE_BEDFILE} > ${OUTPUT_NAME}.allgenenames

cat ${OUTPUT_NAME}.allhconf.genes ${OUTPUT_NAME}.allgenenames | sort | uniq -c | awk '{if ($1==1) {print $2}}' > ${OUTPUT_NAME}.neverfound.genes

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