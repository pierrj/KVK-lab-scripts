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
    bio_rep=$(echo "$line" | cut -f3)
    treatment=$(echo "$line" | cut -f4)
    normalize_factor=$(samtools view -c -F 4 -F 2048 ${normalize_file} | awk -v S=${SCALING_FACTOR} '{print $1/S}' )
    echo -e ${bio_rep} '\t' ${treatment} >> bio_rep_mapfile
    echo ${treatment} >> treatment_mapfile
    awk -v N=${normalize_factor} -v B=${BIN_SIZE} '{sum+=$3} NR%B==0 {print sum/B/N; sum =0}' ${target_file} > ${target_file}.normalized_binned.${bio_rep}.${treatment}
done < ${MAPFILE}

while read line; 
do
    bio_rep=$(echo "$line" | cut -f1)
    treatment=$(echo "$line" | cut -f2)
    paste $(find . -maxdepth 1 -name "*normalized_binned.${bio_rep}*" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${bio_rep}.normalized_binned.${treatment}
done < bio_rep_mapfile

while read treatment; 
do
    paste $(find . -maxdepth 1 -name "*.normalized_binned.${treatment}" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > ${treatment}.normalized_binned
done < treatment_mapfile

## clean up
while read line; 
do
    target_file=$(echo "$line" | cut -f1)
    normalize_file=$(echo "$line" | cut -f2)
    bio_rep=$(echo "$line" | cut -f3)
    treatment=$(echo "$line" | cut -f4)
    mv ${target_file}.normalized_binned.${tech_rep}.${bio_rep}.${treatment} ${target_file}.normalized_binned
done < mapfile

while read line; 
do
    bio_rep=$(echo "$line" | cut -f1)
    treatment=$(echo "$line" | cut -f2)
    mv ${bio_rep}.normalized_binned.${treatment} ${bio_rep}.normalized_binned
done < bio_rep_mapfile

rm bio_rep_mapfile
rm treatment_mapfile