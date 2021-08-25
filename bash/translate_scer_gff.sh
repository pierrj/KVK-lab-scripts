#!/bin/bash


wget http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz

gunzip saccharomyces_cerevisiae.gff.gz

head -28407 saccharomyces_cerevisiae.gff > saccharomyces_cerevisiae.nogenome.gff


if [ -f "scer_translation_column2" ]; then
    rm scer_translation_column2
fi

grep '>' Scer_S288C.fasta | awk -v OFS='\t' '{print substr($1, 2), "chr"substr($7,1,length($7)-1)}' | awk -v OFS='\t' '{if ($2 == "chrgenom") {print $1, "chrmt"} else {print $0}}' > scer_translation


awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{if ($1 !~ /#/) {$1=a[$1];}}1' \
    scer_translation \
    saccharomyces_cerevisiae.nogenome.gff | \
    awk -v OFS='\t' '{ if ($1 !~ /#/) {print $1, $2, $3, $4, $5, $6, $7, $8, $9} else {print $0}}' > Scer_S288C.SGD.translated_to_genbank.gff