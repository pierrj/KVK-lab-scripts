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
while getopts b:e:s:o: option
do
case "${option}"
in
b) BEDFILE=${OPTARG};;
e) ELEMENT=${OPTARG};;
s) SPLIT_READ_FILE=${OPTARG};;
o) OUTPUTNAME=${OPTARG};;
esac
done

## USAGE ##
# this script generate, per ltr retrotransposon in the input, a set of data pertaining to split reads associated with eccdna formation
# -b bedfile of locations for the element in question
# -e name of the element
# -s bamfile containing splitreads
# -o sample name


## move bedfile around here
# filter for specific elements and pieces of elements
grep ${ELEMENT} ${BEDFILE} > tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed
grep LTR tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed > tmp.LTR.${ELEMENT}.${OUTPUTNAME}.filtered.bed
grep INTERNAL tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed > tmp.INTERNAL.${ELEMENT}.${OUTPUTNAME}.filtered.bed

# get list of reads that overlap any ltr portion of the element
## need sort and uniq because some reads overlap two features
bedtools intersect -wa -bed -abam ${SPLIT_READ_FILE} -b tmp.LTR.${ELEMENT}.${OUTPUTNAME}.filtered.bed | sort | uniq > ${ELEMENT}.${OUTPUTNAME}.ltr_overlap_srs.bed
bedtools bamtobed -i ${SPLIT_READ_FILE} > ${OUTPUTNAME}.allsrs.bed

# now that we have the list of reads intersecting ltrs we need to find which secondary location they map to
python /global/home/users/pierrj/git/python/get_other_locs_for_ltr_srs.py ${ELEMENT}.${OUTPUTNAME}.ltr_overlap_srs.bed ${OUTPUTNAME}.allsrs.bed \
    ${ELEMENT}.${OUTPUTNAME}.ltr_and_ltr_overlap_srs.bed ${ELEMENT}.${OUTPUTNAME}.other_overlap_srs.bed

bedtools intersect -wa -a ${ELEMENT}.${OUTPUTNAME}.other_overlap_srs.bed -b tmp.INTERNAL.${ELEMENT}.${OUTPUTNAME}.filtered.bed > ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed

# now grab both locations for each read
awk '{print $4}' ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed > ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.mapfile

if [ -f "${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed" ]; then
    rm ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed
fi

cat ${ELEMENT}.${OUTPUTNAME}.ltr_and_ltr_overlap_srs.bed >> ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed
cat ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed >> ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed

# get just the ltr side of the internal overlap srs
while read sr; do
    grep $sr ${ELEMENT}.${OUTPUTNAME}.ltr_overlap_srs.bed >> ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed
done < ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.mapfile

mv ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed_old
awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6}' ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed_old > ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed

# total counts
ltr_and_ltr_count=$(wc -l ${ELEMENT}.${OUTPUTNAME}.ltr_and_ltr_overlap_srs.bed | awk '{print $1/2}')

ltr_and_internal_count=$(wc -l ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed | awk '{print $1}')

total=$((ltr_and_ltr_count+ltr_and_internal_count))

# then get whole genome coverage for each type of read, split different ways
bedtools coverage -a tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed -b ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed \
 | awk -v OFS='\t' '{print $4, $5}' > ${ELEMENT}.${OUTPUTNAME}.ltr_sr_cov_perfeature

## divided by two bc there are two entries for each split read
bedtools coverage -a tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed -b ${ELEMENT}.${OUTPUTNAME}.ltr_and_ltr_overlap_srs.bed \
 | awk -v OFS='\t' '{print $4, $5/2}' > ${ELEMENT}.${OUTPUTNAME}.ltr_ltr_sr_cov_perfeature

## not divided bc there is only one entry per split read
bedtools coverage -a tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed -b ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed \
 | awk -v OFS='\t' '{print $4, $5}' > ${ELEMENT}.${OUTPUTNAME}.ltr_internal_sr_cov_perfeature

## total read coverage
bedtools coverage -sorted -a tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed \
    -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${OUTPUTNAME}/no_secondary.filtered.sorted.${OUTPUTNAME}.bam | \
    awk -v OFS='\t' '{print $4, $5}' > ${ELEMENT}.${OUTPUTNAME}.read_cov_perfeature

## total number of junction splitreads
bedtools coverage -a tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed \
    -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${OUTPUTNAME}/ltr_eccs_splitreads | \
    awk -v OFS='\t' '{print $4, $5}' > ${ELEMENT}.${OUTPUTNAME}.junction_sr_cov_perfeature

echo -e ${ELEMENT}'\t'${total}