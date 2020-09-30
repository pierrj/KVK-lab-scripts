#!/bin/bash
while getopts p:n:c:s:g: option
do
case "${option}"
in
p) PANNZER_INPUT=${OPTARG};;
n) OUTPUT_NAME=${OPTARG};;
c) SCORE_CUTOFF=${OPTARG};;
m) MAPFILE=${OPTARG};;
b) GENE_BEDFILE=${OPTARG};;
t) TE_BEDFILE=${OPTARG};;
f) PFAM_DIR=${OPTARG};;
a) CDS_FASTA=${OPTARG};;
h) THREADS=${OPTARG};;
esac
done

## rename here

while read SAMPLE
do
    awk '$6=="hconf" && $3-$2>1000' ecccaller_output.${SAMPLE}.renamed.details.tsv > ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv
    line_count=$(wc -l ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv | awk '{print int($1/100)}') 
    sr_count_cutoff=$(sort -k4,4nr ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv | head -${line_count} | tail -1 | awk '{print $4}')
    awk -v CUTOFF=${sr_count_cutoff} '$4>=CUTOFF' ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv > ecccaller_output.${SAMPLE}.common.renamed.details.tsv
done < ${MAPFILE}

## TO GENERATE BIOREP MAPFILE IF STANDARD FILES ARE SET

awk '{print substr($1, 0,4)}' ${MAPFILE} > tmp_biorepmapfile

while read line; 
do
    bio_rep=$(echo "$line" | cut -f1)
    cat $(find . -maxdepth 1 -name "ecccaller_output.${bio_rep}*.common.renamed.details.tsv" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") > ecccaller_output.${biorep}.common.renamed.details.tsv
done < tmp_biorepmapfile

while read line; 
do
    bio_rep=$(echo "$line" | cut -f1)
    bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b ecccaller_output.${bio_rep}.common.renamed.details.tsv | awk '{if ($5!=0) {print $4}}' > ${bio_rep}.common.genes
done < tmp_biorepmapfile

cat $(find . -maxdepth 1 -name "*.common.genes" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | sort | uniq -c | awk '$1==3 {print $2}' > ${OUTPUT_NAME}.common.genes

while read line; 
do
    bio_rep=$(echo "$line" | cut -f1)
    cat $(find . -maxdepth 1 -name "ecccaller_output.${bio_rep}*.hconf.renamed.details.tsv" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") > ${biorep}.hconf.genes
done < tmp_biorepmapfile

cat $(find . -maxdepth 1 -name "${biorep}.hconf.genes" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | sort | uniq -c | awk '$1==3 {print $2}' | sort | uniq > ${OUTPUT_NAME}.hconf.genes

awk '{print $4}' ${GENE_BEDFILE} > ${OUTPUT_NAME}.allgenenames

cat ${OUTPUT_NAME}.hconf ${OUTPUT_NAME}.allgenenames | sort | uniq -c | awk '$1==1' > ${OUTPUT_NAME}.neverfound.genes

Rscript --vanilla /global/home/users/pierrj/git/R/nearest_neighbor_calc_forpipeline.R ${GENE_BEDFILE} ${GENE_BEDFILE} ${OUTPUT_NAME}.common.genes ${OUTPUT_NAME}.common_genes.genedistance.distanceplot.pdf ${OUTPUT_NAME}.common_genes.genedistance.pvalues.txt

Rscript --vanilla /global/home/users/pierrj/git/R/nearest_neighbor_calc_forpipeline.R ${GENE_BEDFILE} ${GENE_BEDFILE} ${OUTPUT_NAME}.common.genes ${OUTPUT_NAME}.neverfound_genes.genedistance.distanceplot.pdf ${OUTPUT_NAME}.neverfound_genes.genedistance.pvalues.txt

Rscript --vanilla /global/home/users/pierrj/git/R/nearest_neighbor_calc_forpipeline.R ${GENE_BEDFILE} ${TE_BEDFILE} ${OUTPUT_NAME}.neverfound.genes ${OUTPUT_NAME}.common_genes.tedistance.distanceplot.pdf ${OUTPUT_NAME}.tedistance.common_genes.pvalues.txt

Rscript --vanilla /global/home/users/pierrj/git/R/nearest_neighbor_calc_forpipeline.R ${GENE_BEDFILE} ${TE_BEDFILE} ${OUTPUT_NAME}.neverfound.genes ${OUTPUT_NAME}.neverfound_genes.tedistance.distanceplot.pdf ${OUTPUT_NAME}.neverfound_genes.common_genes.pvalues.txt

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.common_genes.goenrichment.${SCORE_CUTOFF}_MF -c ${SCORE_CUTOFF} -s ${OUTPUT_NAME}.common.genes -g MF

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.common_genes.goenrichment.${SCORE_CUTOFF}_BP -c ${SCORE_CUTOFF} -s ${OUTPUT_NAME}.common.genes -g BP

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.common_genes.goenrichment.${SCORE_CUTOFF}_CC -c ${SCORE_CUTOFF} -s ${OUTPUT_NAME}.common.genes -g CC

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.neverfound_genes.goenrichment.${SCORE_CUTOFF}_MF -c ${SCORE_CUTOFF} -s ${OUTPUT_NAME}.neverfound.genes -g MF

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.neverfound_genes.goenrichment.${SCORE_CUTOFF}_BP -c ${SCORE_CUTOFF} -s ${OUTPUT_NAME}.neverfound.genes -g BP

/global/home/users/pierrj/git/bash/go_enrichment_wrapper.sh -p ${PANNZER_INPUT} -n ${OUTPUT_NAME}.neverfound_genes.goenrichment.${SCORE_CUTOFF}_CC -c ${SCORE_CUTOFF} -s ${OUTPUT_NAME}.neverfound.genes -g CC

/global/home/users/pierrj/git/bash/pfam_scan_wordcloud_wrapper.sh -g ${OUTPUT_NAME}.common.genes -a ${ALL_GENE_IDS} -n ${OUTPUT_NAME}.common_genes -p ${PFAM_DIR} -c ${CDS_FASTA} -t ${THREADS}

/global/home/users/pierrj/git/bash/pfam_scan_wordcloud_wrapper.sh -g ${OUTPUT_NAME}.neverfound.genes -a ${ALL_GENE_IDS} -n ${OUTPUT_NAME}.never_genes -p ${PFAM_DIR} -c ${CDS_FASTA} -t ${THREADS}