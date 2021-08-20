#!/bin/bash


wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/735/GCA_000001735.2_TAIR10.1/GCA_000001735.2_TAIR10.1_genomic.gff.gz

cat <(awk '{print $1}' /global/scratch/users/pierrj/references/TAIR10_GFF3_genes.gff | sort | uniq | tail -2) <(awk '{print $1}' /global/scratch/users/pierrj/references/TAIR10_GFF3_genes.gff | sort | uniq | head -5)

awk '{ if ($1 !~ /#/) {print $1}}' GCA_000001735.2_TAIR10.1_genomic.gff | sort | uniq > ara_translation_column2

paste ara_translation_column1 ara_translation_column2 > ara_translation

awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{if ($1 !~ /#/) {$1=a[$1];}}1' \
    ara_translation \
    GCA_000001735.2_TAIR10.1_genomic.gff | \
    awk -v OFS='\t' '{ if ($1 !~ /#/) {print $1, $2, $3, $4, $5, $6, $7, $8, $9} else {print $0}}'> GCA_000001735.2_TAIR10.1_genomic.renamed.gff
