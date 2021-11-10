#!/bin/bash
while getopts b:g:o:n: option
do
case "${option}"
in
b) BEDFILE=${OPTARG};;
g) GENOME_FILE=${OPTARG};;
o) OUTPUT_FILE=${OPTARG};;
n) NUM_PRIMER_PAIRS=${OPTARG};;
esac
done

## extend gene entries and get fasta

awk -v OFS='\t' '{print $1, $2-500, $3+500, $4}' $BEDFILE > tmp.extended_bed_coords
bedtools getfasta -fi ${GENOME_FILE} -bed $BEDFILE > tmp.extended_bed_coords.fasta


DOES GET FASTA SPIT OUT SIGNLE LINE OR MULTI LINE FASTA
awk '{ if ($1 !~ />/) {print "SEQUENCE_TEMPLATE="$1}' tmp.extended_bed_coords.fasta > tmp.seq_fastas

awk -v OFS='\t' '{print "SEQUENCE_ID="$4}' ${BEDFILE} > tmp.gene_names

a=($(wc -l tmp.gene_names))

awk -v OFS='\t' '{print "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=0,250,"$3-$2-250",249" }' $BEDFILE > tmp.primer_params

printf "PRIMER_NUM_RETURN=${NUM_PRIMER_PAIRS}\n%.0s" $(eval "echo {1.."$(($a))"}") > tmp.primer_num
printf '=\n%.0s' $(eval "echo {1.."$(($a))"}") > tmp.equal_signs

paste -d"\n" tmp.gene_names tmp.seq_fastas tmp.primer_params tmp.primer_num tmp.equal_signs > tmp.params.boulder
