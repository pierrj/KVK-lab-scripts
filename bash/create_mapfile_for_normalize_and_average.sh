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

## USAGE ##
# this script generates a mapfile to be used with normalize_and_average.sh script
# see description of options above

if [ -f "mapfile_for_normalize_and_average" ]; then
    rm mapfile_for_normalize_and_average
fi

## get sample name from formats

sample=$(head -1 ${MAPFILE} | cut -f1)

start_string_file=$(echo ${FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub(FS $2,"")}1')

end_string_file=$(echo ${FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub($1 FS,"")}1')

if [[ "${NORMALIZE_TYPE}" == "e" ]]
then
# get files names if normalizing to files
start_string_normalize_file=$(echo ${NORMALIZE_FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub(FS $2,"")}1')
end_string_normalize_file=$(echo ${NORMALIZE_FILE} | awk -F"${sample}" 'BEGIN {OFS=FS} {sub($1 FS,"")}1')
elif [[ "${NORMALIZE_TYPE}" != "t" ]]
then
echo "invalid normalize file type"
fi

# write out treatment, biological replicate and technical replicate for each sample in addition to normalization file
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