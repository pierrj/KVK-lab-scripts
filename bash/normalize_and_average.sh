#!/bin/bash
while getopts m:f:b: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
f) SCALING_FACTOR=${OPTARG};; ## for scaling data so numbers arent so gross, usually 1000000
b) BIN_SIZE=${OPTARG};; ## this is for averaging bins to a certain size, should be optional as well, usually 100
esac
done

while read line; 
do
    target_file=$(echo "$line" | cut -f1)
    normalize_file=$(echo "$line" | cut -f2)
    sample=$(echo "$line" | cut -f3)
    bio_rep=$(echo "$line" | cut -f4)
    treatment=$(echo "$line" | cut -f5)
    normalize_factor=$(samtools view -c -F 4 -F 2048 ${normalize_file} | awk -v S=${SCALING_FACTOR} '{print $1/S}' )
    awk -v N=${normalize_factor} -v B=${BIN_SIZE} -v OFS='\t' '{sum+=$3} NR%B==0 {print $1, $2, sum/B/N; sum =0}' ${target_file} > ${sample}.normalized_binned
    cut -f3 ${sample}.normalized_binned > tmp_normalize_and_average.${sample}.normalized_binned.${bio_rep}.${treatment}
done < ${MAPFILE}

awk -v OFS='\t' '{print $1, $2}' ${sample}.normalized_binned > tmp_normalize_and_average_first_two_columns ## THIS LINE MEANS ALL OF THE FILES IN A RUN SHOULD BE MAPPED TO THE SAME GENOME
cut -f4- ${MAPFILE} | sort | uniq > tmp_normalize_and_average_bio_rep_mapfile

while read line; 
do
    bio_rep=$(echo "$line" | cut -f1)
    treatment=$(echo "$line" | cut -f2)
    paste $(find . -maxdepth 1 -name "*normalized_binned.${bio_rep}*" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > tmp_normalize_and_average.${bio_rep}.normalized_binned.${treatment}
    paste tmp_normalize_and_average_first_two_columns tmp_normalize_and_average.${bio_rep}.normalized_binned.${treatment} > ${bio_rep}.normalized_binned 
done < tmp_normalize_and_average_bio_rep_mapfile

cut -f5- ${MAPFILE} | sort | uniq > tmp_normalize_and_average_treatment_mapfile

while read line; 
do
    treatment=$(echo "$line")
    paste $(find . -maxdepth 1 -name "*.normalized_binned.${treatment}" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > tmp_normalize_and_average.${treatment}.normalized_binned
    paste tmp_normalize_and_average_first_two_columns tmp_normalize_and_average.${treatment}.normalized_binned > ${treatment}.normalized_binned 
done < tmp_normalize_and_average_treatment_mapfile

rm tmp_normalize_and_average*