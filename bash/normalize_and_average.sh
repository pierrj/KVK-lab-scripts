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
while getopts m:f:b:c:n: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
f) SCALING_FACTOR=${OPTARG};; ## for scaling data so numbers arent so gross, usually 1000000, if normalize option is n then the input number is not used
b) BIN_SIZE=${OPTARG};; ## this is for averaging bins to a certain size, should be optional as well, usually 100 for bedfiles or 1 for ltr_splitreads
c) COLUMN=${OPTARG};; ## column where the data is (2 for ltr, 3 for genome coverage)
n) NORMALIZE_OPTION=${OPTARG};; ## if f, normalize to file, if n, normalize to number
esac
done

## USAGE ##
# this script is used to average various signals across technical replicates within one biological replicate first and then across all biological replicates
# see description of options above

while read line; 
do
    # grab bio rep and treatment name from input mapfile
    target_file=$(echo "$line" | cut -f1)
    normalize_file=$(echo "$line" | cut -f2)
    sample=$(echo "$line" | cut -f3)
    bio_rep=$(echo "$line" | cut -f4)
    treatment=$(echo "$line" | cut -f5)
    # either table of numbers to normalize too or number of alignments
    if [[ "${NORMALIZE_OPTION}" == "f" ]]
    then
    normalize_factor=$(samtools view -c -F 4 -F 2048 ${normalize_file} | awk -v S=${SCALING_FACTOR} '{print $1/S}' )
    elif [[ "${NORMALIZE_OPTION}" == "n" ]]
    then
    normalize_factor=$(bc -l <<<"${normalize_file}")
    else
    echo "invalid normalize option"
    fi
    # average based off where the number you are averageing is located in the file
    if [[ "${COLUMN}" -eq 3 ]]
    then
    awk -v N=${normalize_factor} -v B=${BIN_SIZE} -v OFS='\t' '{sum+=$3} NR%B==0 {print $1, $2, sum/B/N; sum =0}' ${target_file} > ${sample}.normalized_binned
    cut -f3 ${sample}.normalized_binned > tmp_normalize_and_average.${sample}.normalized_binned.${bio_rep}.${treatment}
    elif [[ "${COLUMN}" -eq 2 ]]
    then
    awk -v N=${normalize_factor} -v B=${BIN_SIZE} -v OFS='\t' '{sum+=$2} NR%B==0 {print $1, sum/B/N; sum =0}' ${target_file} > ${sample}.normalized_binned
    cut -f2 ${sample}.normalized_binned > tmp_normalize_and_average.${sample}.normalized_binned.${bio_rep}.${treatment}
    else
    echo "invalid column number"
    fi
done < ${MAPFILE}

if [[ "${COLUMN}" -eq 3 ]]
then
awk -v OFS='\t' '{print $1, $2}' ${sample}.normalized_binned > tmp_normalize_and_average_first_two_columns ## THIS LINE MEANS ALL OF THE FILES IN A RUN SHOULD BE MAPPED TO THE SAME GENOME
cut -f4- ${MAPFILE} | sort | uniq > tmp_normalize_and_average_bio_rep_mapfile
elif [[ "${COLUMN}" -eq 2 ]]
then
awk -v OFS='\t' '{print $1}' ${sample}.normalized_binned > tmp_normalize_and_average_first_column ## THIS LINE MEANS ALL OF THE FILES IN A RUN SHOULD BE MAPPED TO THE SAME GENOME
cut -f4- ${MAPFILE} | sort | uniq > tmp_normalize_and_average_bio_rep_mapfile
else
echo "invalid column number"
fi

## average across biological replicates based off the names of the files
if [[ "${COLUMN}" -eq 3 ]]
then
while read line; 
do
    bio_rep=$(echo "$line" | cut -f1)
    treatment=$(echo "$line" | cut -f2)
    paste $(find . -maxdepth 1 -name "*.normalized_binned.${bio_rep}*" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > tmp_normalize_and_average.average.${bio_rep}.normalized_binned.${treatment}
    # calculate standard deviation as well
    paste $(find . -maxdepth 1 -name "*.normalized_binned.${bio_rep}*" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{ A=0; V=0; for(N=1; N<=NF; N++) A+=$N ; A/=NF ; for(N=1; N<=NF; N++) V+=(($N-A)*($N-A))/(NF-1); print sqrt(V) }' > tmp_normalize_and_average.stdev.${bio_rep}.normalized_binned.${treatment}
    paste tmp_normalize_and_average_first_two_columns tmp_normalize_and_average.average.${bio_rep}.normalized_binned.${treatment} tmp_normalize_and_average.stdev.${bio_rep}.normalized_binned.${treatment} > ${bio_rep}.normalized_binned 
done < tmp_normalize_and_average_bio_rep_mapfile
elif [[ "${COLUMN}" -eq 2 ]]
then
while read line; 
do
    bio_rep=$(echo "$line" | cut -f1)
    treatment=$(echo "$line" | cut -f2)
    paste $(find . -maxdepth 1 -name "*.normalized_binned.${bio_rep}*" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > tmp_normalize_and_average.average.${bio_rep}.normalized_binned.${treatment}
    # calculate standard deviation as well
    paste $(find . -maxdepth 1 -name "*.normalized_binned.${bio_rep}*" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{ A=0; V=0; for(N=1; N<=NF; N++) A+=$N ; A/=NF ; for(N=1; N<=NF; N++) V+=(($N-A)*($N-A))/(NF-1); print sqrt(V) }' > tmp_normalize_and_average.stdev.${bio_rep}.normalized_binned.${treatment}
    paste tmp_normalize_and_average_first_column tmp_normalize_and_average.average.${bio_rep}.normalized_binned.${treatment} tmp_normalize_and_average.stdev.${bio_rep}.normalized_binned.${treatment} > ${bio_rep}.normalized_binned 
done < tmp_normalize_and_average_bio_rep_mapfile
else
echo "invalid column number"
fi

# make mapfile for second round of averaging
cut -f5- ${MAPFILE} | sort | uniq > tmp_normalize_and_average_treatment_mapfile

if [[ "${COLUMN}" -eq 3 ]]
then
while read line; 
do
    treatment=$(echo "$line")
    paste $(find . -maxdepth 1 -name "*.average.*.normalized_binned.${treatment}" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > tmp_normalize_and_average.average.${treatment}.normalized_binned
    # calculate standard deviation as well    
    paste $(find . -maxdepth 1 -name "*.stdev.*.normalized_binned.${treatment}" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i^2; sum /= NF; print sqrt(sum)}' > tmp_normalize_and_average.stdev.${treatment}.normalized_binned
    paste tmp_normalize_and_average_first_two_columns tmp_normalize_and_average.average.${treatment}.normalized_binned tmp_normalize_and_average.stdev.${treatment}.normalized_binned > ${treatment}.normalized_binned 
done < tmp_normalize_and_average_treatment_mapfile
elif [[ "${COLUMN}" -eq 2 ]]
then
while read line; 
do
    treatment=$(echo "$line")
    paste $(find . -maxdepth 1 -name "*.average.*.normalized_binned.${treatment}" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > tmp_normalize_and_average.average.${treatment}.normalized_binned
    # calculate standard deviation as well    
    paste $(find . -maxdepth 1 -name "*.stdev.*.normalized_binned.${treatment}" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i^2; sum /= NF; print sqrt(sum)}' > tmp_normalize_and_average.stdev.${treatment}.normalized_binned
    paste tmp_normalize_and_average_first_column tmp_normalize_and_average.average.${treatment}.normalized_binned tmp_normalize_and_average.stdev.${treatment}.normalized_binned > ${treatment}.normalized_binned 
done < tmp_normalize_and_average_treatment_mapfile
else
echo "invalid column number"
fi

rm tmp_normalize_and_average*