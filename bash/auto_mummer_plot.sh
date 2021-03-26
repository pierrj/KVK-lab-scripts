#!/bin/bash
while getopts r:g:s:p: option
do
case "${option}"
in
r) REFERENCE=${OPTARG};;
g) GENOME_DIR=${OPTARG};;
s) SUBDIR_OUTPUT=${OPTARG};;
p) PERCENT_ZEROES_FILTER=${OPTARG};;
esac
done

if [ -d "${SUBDIR_OUTPUT}" ]; then
    rm -r ${SUBDIR_OUTPUT}
fi

ls -1 | grep -v genome_mapfile | awk '{print substr($1,9)}'> genome_mapfile

mkdir ${SUBDIR_OUTPUT}

while read genome ; do
    cd guy11_v_${genome}
    SUBDIR_TMP=tmp
    if [ -d "${SUBDIR_TMP}" ]; then
        rm -r ${SUBDIR_TMP}
    fi
    mkdir ${SUBDIR_TMP}
    alignment_number=$(ls -1 | wc -l | awk '{print ($1-1)/2-1}')
    seq 0 ${alignment_number} > ${SUBDIR_TMP}/mapfile
    ISOLATE=${GENOME_DIR}/${genome}_genomic.fna
    while read plot; do
        bedtools getfasta -fi ${REFERENCE} -bed guy11_${plot}.bed -fo ${SUBDIR_TMP}/guy11_${plot}.fasta
        samtools faidx ${SUBDIR_TMP}/guy11_${plot}.fasta
        cut -f1,2 ${SUBDIR_TMP}/guy11_${plot}.fasta.fai > ${SUBDIR_TMP}/guy11_${plot}.genomesize
        bedtools getfasta -fi ${ISOLATE} -bed ${genome}_${plot}.bed -fo ${SUBDIR_TMP}/${genome}_${plot}.fasta
        nucmer -p ${SUBDIR_TMP}/${plot} ${SUBDIR_TMP}/guy11_${plot}.fasta ${SUBDIR_TMP}/${genome}_${plot}.fasta
        show-coords ${SUBDIR_TMP}/${plot}.delta | tail -n +6 | awk -v OFS='\t' '{print $12, $1, $2}' | sort -k1,1 -k2,2n >  ${SUBDIR_TMP}/${plot}.bed
        bedtools genomecov -d -i ${SUBDIR_TMP}/${plot}.bed -g ${SUBDIR_TMP}/guy11_${plot}.genomesize > ${SUBDIR_TMP}/${plot}.genomecov
        total_size=$(wc -l ${SUBDIR_TMP}/${plot}.genomecov | awk '{print $1}')
        size_zeroes=$(awk '$3==0' ${SUBDIR_TMP}/${plot}.genomecov | wc -l | awk '{print $1}')
        percent_zeroes=$(awk -v var1=$size_zeroes -v var2=$total_size 'BEGIN { OFMT="%f";print  ( var1 / var2 ) }')
        if (( $(echo "$percent_zeroes < ${PERCENT_ZEROES_FILTER}" |bc -l) ))
        then
            mummerplot --color -postscript -p ${SUBDIR_TMP}/${plot} ${SUBDIR_TMP}/${plot}.delta
            ps2pdf ${SUBDIR_TMP}/${plot}.ps ${SUBDIR_TMP}/${plot}.pdf
            convert -density 150 ${SUBDIR_TMP}/${plot}.pdf -quality 90 ../${SUBDIR_OUTPUT}/${genome}_${plot}.jpg
        fi
    done < ${SUBDIR_TMP}/mapfile
    cd ..
done < genome_mapfile