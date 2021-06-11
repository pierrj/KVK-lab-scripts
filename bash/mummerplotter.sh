#!/bin/bash
while getopts r:q:e:u:o:p: option
do
case "${option}"
in
r) REFERENCE=${OPTARG};;
q) QUERY=${OPTARG};;
e) ref=${OPTARG};;
u) quer=${OPTARG};;
o) SUBDIR_OUTPUT=${OPTARG};;
p) PERCENT_ZEROES_FILTER=${OPTARG};;
esac
done

output_realpath=$(realpath ${SUBDIR_OUTPUT})
SUBDIR_TMP=tmp
if [ -d "${SUBDIR_TMP}" ]; then
    rm -r ${SUBDIR_TMP}
fi
mkdir ${SUBDIR_TMP}
alignment_number=$(ls -1 | wc -l | awk '{print ($1-1)/2-1}')
seq 0 ${alignment_number} > ${SUBDIR_TMP}/mapfile
while read plot; do
    bedtools getfasta -fi ${REFERENCE} -bed ${ref}_${plot}.bed -fo ${SUBDIR_TMP}/${ref}_${plot}.fasta
    samtools faidx ${SUBDIR_TMP}/${ref}_${plot}.fasta
    cut -f1,2 ${SUBDIR_TMP}/${ref}_${plot}.fasta.fai > ${SUBDIR_TMP}/${ref}_${plot}.genomesize
    bedtools getfasta -fi ${QUERY} -bed ${quer}_${plot}.bed -fo ${SUBDIR_TMP}/${quer}_${plot}.fasta
    /global/scratch/pierrj/mummer_4/bin/nucmer -p ${SUBDIR_TMP}/${plot} ${SUBDIR_TMP}/${ref}_${plot}.fasta ${SUBDIR_TMP}/${quer}_${plot}.fasta
    /global/scratch/pierrj/mummer_4/bin/show-coords ${SUBDIR_TMP}/${plot}.delta | tail -n +6 | awk -v OFS='\t' '{print $12, $1, $2}' | sort -k1,1 -k2,2n >  ${SUBDIR_TMP}/${plot}.bed
    bedtools genomecov -d -i ${SUBDIR_TMP}/${plot}.bed -g ${SUBDIR_TMP}/${ref}_${plot}.genomesize > ${SUBDIR_TMP}/${plot}.genomecov
    total_size=$(wc -l ${SUBDIR_TMP}/${plot}.genomecov | awk '{print $1}')
    size_zeroes=$(awk '$3==0' ${SUBDIR_TMP}/${plot}.genomecov | wc -l | awk '{print $1}')
    percent_zeroes=$(awk -v var1=$size_zeroes -v var2=$total_size 'BEGIN { OFMT="%f";print  ( var1 / var2 ) }')
    if (( $(echo "$percent_zeroes < ${PERCENT_ZEROES_FILTER}" |bc -l) ))
    then
        /global/scratch/pierrj/mummer_4/bin/mummerplot --color -postscript -p ${SUBDIR_TMP}/${plot} ${SUBDIR_TMP}/${plot}.delta
        /global/home/users/pierrj/ps2pdf/usr/bin/ps2pdf ${SUBDIR_TMP}/${plot}.ps ${SUBDIR_TMP}/${plot}.pdf
        convert -density 150 ${SUBDIR_TMP}/${plot}.pdf -quality 90 ${output_realpath}/${ref}_v_${quer}_${plot}.jpg
        cp ${SUBDIR_TMP}/${plot}.pdf ${output_realpath}/${ref}_v_${quer}_${plot}.pdf
    fi
done < ${SUBDIR_TMP}/mapfile