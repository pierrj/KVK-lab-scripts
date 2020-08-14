#!/bin/bash
## input file name/options
while read SAMPLE
do
    echo $SAMPLE
    wc -l ecccaller_output.${SAMPLE}.renamed.details.tsv
    awk '$6=="hconf" && $3-$2>1000' ecccaller_output.${SAMPLE}.renamed.details.tsv > ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv
    wc -l ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv
    line_count=$(wc -l ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv | awk '{print int($1/100)}') 
    echo $line_count
    sr_count_cutoff=$(sort -k4,4nr ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv | head -${line_count} | tail -1 | awk '{print $4}')
    echo $sr_count_cutoff
    awk -v CUTOFF=${sr_count_cutoff} '$4>=CUTOFF' ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv > ecccaller_output.${SAMPLE}.common.renamed.details.tsv
    wc -l ecccaller_output.${SAMPLE}.common.renamed.details.tsv
done < mapfile

cat ecccaller_output.G3_1A.common.renamed.details.tsv ecccaller_output.G3_1B.common.renamed.details.tsv ecccaller_output.G3_1C.common.renamed.details.tsv > ecccaller_output.G3_1.common.renamed.details.tsv

cat ecccaller_output.G3_2A.common.renamed.details.tsv ecccaller_output.G3_2B.common.renamed.details.tsv ecccaller_output.G3_2C.common.renamed.details.tsv > ecccaller_output.G3_2.common.renamed.details.tsv

cat ecccaller_output.G3_3A.common.renamed.details.tsv ecccaller_output.G3_3B.common.renamed.details.tsv > ecccaller_output.G3_3.common.renamed.details.tsv

cd genes_on_eccs

bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.just_genes.bed -b ../ecccaller_output.G3_1.common.renamed.details.tsv | awk '{if ($5!=0) {print $4}}' > G3_1.common.genes

bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.just_genes.bed -b ../ecccaller_output.G3_2.common.renamed.details.tsv | awk '{if ($5!=0) {print $4}}' > G3_2.common.genes

bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.just_genes.bed -b ../ecccaller_output.G3_3.common.renamed.details.tsv | awk '{if ($5!=0) {print $4}}' > G3_3.common.genes

wc -l G3_1.common.genes

wc -l G3_2.common.genes

wc -l G3_3.common.genes

cat G3_1.common.genes G3_3.common.genes G3_2.common.genes | sort | uniq -c | awk '$1==3 {print $2}' > G3_all.common.genes