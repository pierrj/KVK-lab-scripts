#!/bin/bash
while getopts s:b:p:g: option
do
case "${option}"
in
s) SAMPLE=${OPTARG};;
b) FILTERED_BAMFILE=${OPTARG};;## THIS SHOULD BE COORDINATE SORTED AND INCLUDE ONLY ONE RECORD PER ALIGNMENT AND IT SHOULD NOT INCLUDE MAPQ0 READS
p) PACBIO_MAPPED=${OPTARG};;
g) GENE_BEDFILE=${OPTARG};;
esac
done

samtools view -F 4 ${PACBIO_MAPPED} | awk '{print $2}' > ${SAMPLE}_aligned_samflags
bedtools bamtobed -cigar -i ${PACBIO_MAPPED} | paste - ${SAMPLE}_aligned_samflags > ${SAMPLE}_aligned_pacbio.bed

samtools view -f 65 -F 4 ${FILTERED_BAMFILE} > tmp.1.${SAMPLE}.sam
splitread_file="1.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} > tmp.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{ if (and($2, 16) && $18 == "1")
        {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, 2}
    else if (and($2, 16) && $18 == "2")
        {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, 1}
    else
        {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18}
}' tmp.exactlytwice.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.startendsorted.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f3=$3
    getline 
    if (f3 != $3) {
        print prev
        print $0
    }
}' tmp.startendsorted.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.diff_chrom.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 129 -F 4 ${FILTERED_BAMFILE} > tmp.2.${SAMPLE}.sam
splitread_file="2.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} > tmp.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{ if (and($2, 16) && $18 == "1")
        {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, 2}
    else if (and($2, 16) && $18 == "2")
        {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, 1}
    else
        {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18}
}' tmp.exactlytwice.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.startendsorted.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f3=$3
    getline 
    if (f3 != $3) {
        print prev
        print $0
    }
}' tmp.startendsorted.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.diff_chrom.exactlytwice.qualityfiltered.${splitread_file}

cat tmp.diff_chrom.exactlytwice.qualityfiltered.1.${SAMPLE}.sam tmp.diff_chrom.exactlytwice.qualityfiltered.2.${SAMPLE}.sam > tmp.diff_chrom.exactlytwice.qualityfiltered.${SAMPLE}.sam
python /global/home/users/pierrj/git/python/filter_for_match_lengths_adjustable.py tmp.diff_chrom.exactlytwice.qualityfiltered.${SAMPLE}.sam tmp.match_length_filtered.diff_chrom.exactlytwice.qualityfiltered.${SAMPLE}.sam 20
samtools view -b -h <(cat <(samtools view -H ${FILTERED_BAMFILE}) tmp.match_length_filtered.diff_chrom.exactlytwice.qualityfiltered.${SAMPLE}.sam) > tmp.match_length_filtered.diff_chrom.exactlytwice.qualityfiltered.${SAMPLE}.bam
bedtools bamtobed -i tmp.match_length_filtered.diff_chrom.exactlytwice.qualityfiltered.${SAMPLE}.bam | awk -v OFS='\t' '{
    prev=$0; f1=$1; f2=$2; f3=$3; f4=$4; f6=$6
    getline 
    if ($4 == f4 && f6 == "+" && $6 == "+") {
        print f1, f3, $1, $2, $4
    }
    else if ($4 == f4 && f6 == "+" && $6 == "-") {
        print f1, f3, $1, $3, $4
    }
    else if ($4 == f4 && f6 == "-" && $6 == "+") {
        print f1, f2, $1, $2, $4
    }
    else if ($4 == f4 && f6 == "-" && $6 == "-") {
        print f1, f2, $1, $3, $4
    }
}' > ${SAMPLE}.chimeric_ecc_splitreads.bed

python /global/home/users/pierrj/git/python/get_chimeric_eccdnas_with_genes.py ${SAMPLE}_aligned_pacbio.bed ${SAMPLE}.chimeric_ecc_splitreads.bed ${SAMPLE}_chimeric_eccs_partial_circles ${SAMPLE}_chimeric_eccs_full_circles ${SAMPLE}_chimeric_junctions_with_genes ${SAMPLE} 20 50000 ${GENE_BEDFILE}