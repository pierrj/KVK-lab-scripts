#!/bin/bash
#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
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

## USAGE ##
# this script takes a bed file of eccdna locations and designs primers to validate those junctions
# -b bedfile of eccdna locations
# -g genome fasta file
# -o output file for primer pairs
# -n number of primer pairs to generate

## get junction fasta and turn into primer3 input boulder file
## bedfile should have circle names in fourth column
awk -v OFS='\t' '{print $1, $2, $2+200, $4; print $1, $3-200, $3, $4}' ${BEDFILE} > tmp.junction_coords.bed
bedtools getfasta -fi ${GENOME_FILE} -bed tmp.junction_coords.bed > tmp.junctions.fasta
while read -r ONE; do read -r TWO; read -r THREE; read -r FOUR; echo "SEQUENCE_TEMPLATE=${FOUR}${TWO}"; done < tmp.junctions.fasta > tmp.merged_fastas
awk -v OFS='\t' '{print "SEQUENCE_ID="$4}' ${BEDFILE} > tmp.circle_names
a=($(wc -l tmp.circle_names))
## primer3 params
printf 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=0,195,205,195\n%.0s' $(eval "echo {1.."$(($a))"}") > tmp.primer_params
printf "PRIMER_NUM_RETURN=${NUM_PRIMER_PAIRS}\n%.0s" $(eval "echo {1.."$(($a))"}") > tmp.primer_num
printf '=\n%.0s' $(eval "echo {1.."$(($a))"}") > tmp.equal_signs
paste -d"\n" tmp.circle_names tmp.merged_fastas tmp.primer_params tmp.primer_num tmp.equal_signs > tmp.junctions.boulder

# generate primers through primer 3
/global/scratch/users/pierrj/scripts/primer3/src/primer3_core --p3_settings_file=/global/scratch/users/pierrj/scripts/primer3/src/PMJ_settings_junction_primers tmp.junctions.boulder > tmp.primer3.boulder
grep -oP '(PRIMER_.*._.*._SEQUENCE=)\K.*' tmp.primer3.boulder > tmp.junctions.primer_seqs
awk '{print $4}' ${BEDFILE} > tmp.circle_names_awk

if [ -f "tmp.primer_names" ]; then
    rm tmp.primer_names
fi

# grab primer names
while read circle; do
    for i in $(seq 1 $NUM_PRIMER_PAIRS); do
        echo ${circle}_F${i} >> tmp.primer_names
        echo ${circle}_R${i} >> tmp.primer_names
    done
done < tmp.circle_names_awk

if [ -f "tmp.primer_names_product_size" ]; then
    rm tmp.primer_names_product_size
fi

# grab product sizes
while read circle; do
    for i in $(seq 1 $NUM_PRIMER_PAIRS); do
        echo ${circle}_F${i}_${circle}_R${i} >> tmp.primer_names_product_size
    done
done < tmp.circle_names_awk

## make output file for export into IDT for ordering
paste -d';' tmp.primer_names tmp.junctions.primer_seqs > ${OUTPUT_FILE}
## also output file of output sizes
grep -oP '(PRIMER_PAIR_.*._PRODUCT_SIZE=)\K.*' tmp.primer3.boulder | paste tmp.primer_names_product_size - > ${OUTPUT_FILE}.product_sizes