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


## move bedfile around here
grep ${ELEMENT} ${BEDFILE} > tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed
grep LTR tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed > tmp.LTR.${ELEMENT}.${OUTPUTNAME}.filtered.bed
grep INTERNAL tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed > tmp.INTERNAL.${ELEMENT}.${OUTPUTNAME}.filtered.bed

## need sort and uniq because some reads overlap two features
## maybe this means I need to differentiate between actual split reads and reads overlapping two features right next to each other?
bedtools intersect -wa -bed -abam ${SPLIT_READ_FILE} -b tmp.LTR.${ELEMENT}.${OUTPUTNAME}.filtered.bed | sort | uniq > ${ELEMENT}.${OUTPUTNAME}.ltr_overlap_srs.bed

bedtools bamtobed -i ${SPLIT_READ_FILE} > ${OUTPUTNAME}.allsrs.bed

python /global/home/users/pierrj/git/python/get_other_locs_for_ltr_srs.py ${ELEMENT}.${OUTPUTNAME}.ltr_overlap_srs.bed ${OUTPUTNAME}.allsrs.bed ${ELEMENT}.${OUTPUTNAME}.ltr_and_ltr_overlap_srs.bed ${ELEMENT}.${OUTPUTNAME}.other_overlap_srs.bed

bedtools intersect -wa -a ${ELEMENT}.${OUTPUTNAME}.other_overlap_srs.bed -b tmp.INTERNAL.${ELEMENT}.${OUTPUTNAME}.filtered.bed > ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed

ltr_and_ltr_count=$(wc -l ${ELEMENT}.${OUTPUTNAME}.ltr_and_ltr_overlap_srs.bed | awk '{print $1/2}')

ltr_and_internal_count=$(wc -l ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed | awk '{print $1}')

total=$((ltr_and_ltr_count+ltr_and_internal_count))

echo -e ${ELEMENT}'\t'${total}