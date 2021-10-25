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

# generate primers
/global/scratch/users/pierrj/scripts/primer3/src/primer3_core --p3_settings_file=/global/scratch/users/pierrj/scripts/primer3/src/PMJ_settings_junction_primers tmp.junctions.boulder > tmp.primer3.boulder
grep -oP '(PRIMER_.*._.*._SEQUENCE=)\K.*' tmp.primer3.boulder > tmp.junctions.primer_seqs
awk '{print $4}' ${BEDFILE} > tmp.circle_names_awk

if [ -f "tmp.primer_names" ]; then
    rm tmp.primer_names
fi

while read circle; do
    for i in $(seq 1 $NUM_PRIMER_PAIRS); do
        echo ${circle}_F${i} >> tmp.primer_names
        echo ${circle}_R${i} >> tmp.primer_names
    done
done < tmp.circle_names_awk

if [ -f "tmp.primer_names_product_size" ]; then
    rm tmp.primer_names_product_size
fi

while read circle; do
    for i in $(seq 1 $NUM_PRIMER_PAIRS); do
        echo ${circle}_F${i}_${circle}_R${i} >> tmp.primer_names_product_size
    done
done < tmp.circle_names_awk

## make output file for export into IDT
paste -d';' tmp.primer_names tmp.junctions.primer_seqs > ${OUTPUT_FILE}
grep -oP '(PRIMER_PAIR_.*._PRODUCT_SIZE=)\K.*' tmp.primer3.boulder | paste tmp.primer_names_product_size - > ${OUTPUT_FILE}.product_sizes