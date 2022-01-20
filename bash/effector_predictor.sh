#!/bin/bash
while getopts i:o: option
do
case "${option}"
in
i) INPUT_FILE=${OPTARG};;
o) OUTPUT_NAME=${OPTARG};;
esac
done

cat ${INPUT_FILE} | /global/scratch/users/pierrj/tmhmm-2.0c/bin/tmhmm > ${OUTPUT_NAME}.tmhmmout

awk '{if ($7 == 0){print $2}} ' ${OUTPUT_NAME}.tmhmmout > ${OUTPUT_NAME}.notm.names

seqtk subseq ${INPUT_FILE} ${OUTPUT_NAME}.notm.names > ${OUTPUT_NAME}.notm.faa

/global/scratch/users/pierrj/signalp-4.1/signalp -t euk -u 0.34 -U 0.34 -f short ${OUTPUT_NAME}.notm.faa > ${OUTPUT_NAME}.signalpout.short

awk '{ if ($10 == "Y") {print $1} }' ${OUTPUT_NAME}.signalpout.short > ${OUTPUT_NAME}.sp.names

seqtk subseq ${OUTPUT_NAME}.notm.faa ${OUTPUT_NAME}.sp.names > ${OUTPUT_NAME}.notm.sp.faa

python /global/scratch/users/pierrj/EffectorP-3.0/EffectorP.py -i ${OUTPUT_NAME}.notm.sp.faa > ${OUTPUT_NAME}.effectorpout

awk '{ if ($(NF) == "effector") {print $1}}' ${OUTPUT_NAME}.effectorpout > ${OUTPUT_NAME}.predicted_effectors