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
### this is fine i think
bedtools intersect -wa -bed -abam ${SPLIT_READ_FILE} -b tmp.LTR.${ELEMENT}.${OUTPUTNAME}.filtered.bed | sort | uniq > ${ELEMENT}.${OUTPUTNAME}.ltr_overlap_srs.bed

## MAKE SURE THE FACT THAT READS 1 AND READS2 ARENT MERGED TOGETHER AT SOME POINT HERE!!!!
bedtools bamtobed -i ${SPLIT_READ_FILE} > ${OUTPUTNAME}.allsrs.bed

python /global/home/users/pierrj/git/python/get_other_locs_for_ltr_srs.py ${ELEMENT}.${OUTPUTNAME}.ltr_overlap_srs.bed ${OUTPUTNAME}.allsrs.bed \
    ${ELEMENT}.${OUTPUTNAME}.ltr_and_ltr_overlap_srs.bed ${ELEMENT}.${OUTPUTNAME}.other_overlap_srs.bed

bedtools intersect -wa -a ${ELEMENT}.${OUTPUTNAME}.other_overlap_srs.bed -b tmp.INTERNAL.${ELEMENT}.${OUTPUTNAME}.filtered.bed > ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed

awk '{print $4}' ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed > ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.mapfile

if [ -f "${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed" ]; then
    rm ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed
fi

cat ${ELEMENT}.${OUTPUTNAME}.ltr_and_ltr_overlap_srs.bed >> ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed
cat ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed >> ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed

while read sr; do
    grep $sr ${ELEMENT}.${OUTPUTNAME}.ltr_overlap_srs.bed >> ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed
done < ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.mapfile

mv ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed_old
awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6}' ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed_old > ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed

ltr_and_ltr_count=$(wc -l ${ELEMENT}.${OUTPUTNAME}.ltr_and_ltr_overlap_srs.bed | awk '{print $1/2}')

ltr_and_internal_count=$(wc -l ${ELEMENT}.${OUTPUTNAME}.ltr_and_internal_overlap_srs.bed | awk '{print $1}')

total=$((ltr_and_ltr_count+ltr_and_internal_count))

bedtools coverage -a tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed -b ${ELEMENT}.${OUTPUTNAME}.all_ltr_sr_locations.bed \
 | awk -v OFS='\t' '{print $4, $5}' > ${ELEMENT}.${OUTPUTNAME}.ltr_sr_cov_perfeature

bedtools coverage -sorted -a tmp.${ELEMENT}.${OUTPUTNAME}.filtered.bed \
    -b /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${OUTPUTNAME}/no_secondary.filtered.sorted.${OUTPUTNAME}.bam > ${ELEMENT}.${OUTPUTNAME}.read_cov_perfeature

echo -e ${ELEMENT}'\t'${total}