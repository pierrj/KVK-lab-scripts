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
while getopts b:m:s:o: option
do
case "${option}"
in
b) BEDFILE=${OPTARG};;
m) MAPFILE=${OPTARG};;
s) SPLIT_READ_FILE=${OPTARG};;
o) OUTPUTNAME=${OPTARG};;
esac
done

## USAGE ##
# this script takes a bam file containing split reads of interest that might be related to ltr retrotransposon and associates them to a list of ltr retrotransposons
# options:
# -b bedfile of te sequences
# -m mapfile containing te names
# -s split read bam file
# -o name of output

if [ -f "${OUTPUTNAME}.sr_count_per_element" ]; then
    rm ${OUTPUTNAME}.sr_count_per_element
fi

if [ -f "${OUTPUTNAME}.ltr_sr_cov_perfeature" ]; then
    rm ${OUTPUTNAME}.ltr_sr_cov_perfeature
fi

if [ -f "${OUTPUTNAME}.read_cov_perfeature" ]; then
    rm ${OUTPUTNAME}.read_cov_perfeature
fi


# loop through elements and get different types of split reads and read coverage information per transposable element
# see get_ltr_sr_circle_count_per_element.sh for more detail
while read element
do
    /global/home/users/pierrj/git/bash/get_ltr_sr_circle_count_per_element.sh -b ${BEDFILE} -e ${element} -s ${SPLIT_READ_FILE} -o ${OUTPUTNAME} >> ${OUTPUTNAME}.sr_count_per_element
    cat ${element}.${OUTPUTNAME}.ltr_sr_cov_perfeature >> ${OUTPUTNAME}.ltr_sr_cov_perfeature
    cat ${element}.${OUTPUTNAME}.ltr_ltr_sr_cov_perfeature >> ${OUTPUTNAME}.ltr_ltr_sr_cov_perfeature
    cat ${element}.${OUTPUTNAME}.ltr_internal_sr_cov_perfeature >> ${OUTPUTNAME}.ltr_internal_sr_cov_perfeature
    cat ${element}.${OUTPUTNAME}.read_cov_perfeature >> ${OUTPUTNAME}.read_cov_perfeature
    cat ${element}.${OUTPUTNAME}.junction_sr_cov_perfeature >> ${OUTPUTNAME}.junction_sr_cov_perfeature
done < ${MAPFILE}