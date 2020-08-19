#!/bin/bash
while getopts m:s:t:b: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
b) FILTERED_BAMFILE=${OPTARG};;
esac
done

if [ -f "${SAMPLE}.rarefaction_output_hconf" ]; then
    rm ${SAMPLE}.rarefaction_output_hconf
fi

if [ -f "${SAMPLE}.rarefaction_output_all" ]; then
    rm ${SAMPLE}.rarefaction_output_all
fi

for i in $(seq 0.1 0.1 1.0); do
    mkdir $i
    cd $i
    rarefaction_sample=tmp.rarefaction.${SAMPLE}_${i}
    samtools view -bh -s ${i} ${FILTERED_BAMFILE} > ${rarefaction_sample}.bam
    /global/home/users/pierrj/git/bash/call_ecc_regions.sh -m ${MAPFILE} -s ${rarefaction_sample} -t ${THREADS} -b ${rarefaction_sample}.bam
    /global/home/users/pierrj/git/bash/assign_confidence_nodb.sh -m ${MAPFILE} -s ${rarefaction_sample} -t ${THREADS} -b ${rarefaction_sample}.bam -r ${rarefaction_sample}.confirmedsplitreads.bed
    awk '{ if ($6 == "hconf") print $0}' ecccaller_output.${rarefaction_sample}.details.tsv | wc -l | awk '{print $1}' >> ../${SAMPLE}.rarefaction_output_hconf
    wc -l ecccaller_output.${rarefaction_sample}.details.tsv | awk '{print $1}' >> ../${SAMPLE}.rarefaction_output_all
    cd ..
done

# rm *tmp.rarefaction.*