#!/bin/bash
while getopts m:t:n:y: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
t) FILE=${OPTARG};; # example file format
n) NORMALIZE_FILE=${OPTARG};; ## example file if y option is e, otherwise table of numbers if option is t
y) NORMALIZE_TYPE=${OPTARG};; ## can be either an example file name, e, or a table of numbers, t
esac
done

if [ -f "mapfile_for_normalize_and_average" ]; then
    rm mapfile_for_normalize_and_average
fi

## GET SAMPLE ##

sample=$(head -1 ${MAPFILE} | cut -f1)

start_string_file=$(echo ${FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub(FS $2,"")}1')

end_string_file=$(echo ${FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub($1 FS,"")}1')

if [[ "${NORMALIZE_TYPE}" == "e" ]]
then
start_string_normalize_file=$(echo ${NORMALIZE_FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub(FS $2,"")}1')
end_string_normalize_file=$(echo ${NORMALIZE_FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub($1 FS,"")}1')
elif [[ "${NORMALIZE_TYPE}" != "t" ]]
then
echo "invalid normalize file type"
fi

while read line; 
do
    sample=$(echo "$line" | cut -f1)
    bio_rep=$(echo "$line" | cut -f2)
    treatment=$(echo "$line" | cut -f3)
    cd ${sample}
    file_path=$(realpath ${start_string_file}${sample}${end_string_file})
    if [[ "${NORMALIZE_TYPE}" == "e" ]]
    then
    normalize_path=$(realpath ${start_string_normalize_file}${sample}${end_string_normalize_file})
    fi
    cd ..
    if [[ "${NORMALIZE_TYPE}" == "t" ]]
    then
    normalize_path=$(grep ${sample} ${NORMALIZE_FILE} | awk '{print $2}')
    fi
    echo -e ${file_path}'\t' ${normalize_path}'\t'${sample}'\t'${bio_rep}'\t'${treatment} >> mapfile_for_normalize_and_average
done < ${MAPFILE}