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
while getopts m: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
esac
done

## USAGE ##
# this script inputs a list of sample names and breaks them down into a table of samples and the treatments, biorep and tech reps that are associated with them

if [ -f "sample_mapfile" ]; then
    rm sample_mapfile
fi

while read sample; 
do
    bio_rep=$(echo ${sample} | cut -c-4)
    treatment=$(echo ${sample} | cut -c-2)
    echo -e ${sample}'\t'${bio_rep}'\t'${treatment} >> sample_mapfile
done < ${MAPFILE}