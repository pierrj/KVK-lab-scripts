#!/bin/bash
while getopts b:d:g:w:t:o:s: option
do
case "${option}"
in
b) REGIONS_BED=${OPTARG};;
d) DENSITY_FILE+=(${OPTARG});; ## either one bed file or two bam files
g) GENOME_FILE=${OPTARG};;
w) WINDOWS=${OPTARG};;
t) THREADS=${OPTARG};;
o) OUTPUT_NAME=${OPTARG};;
s) SV=${OPTARG};;
esac
done

genome_basename=$(basename ${GENOME_FILE})
genome_basename=${genome_basename%%.*}

if [ ! -f "${GENOME_FILE}.chromsizes" ]; then
    samtools faidx ${GENOME_FILE}
    cut -f1,2 ${GENOME_FILE}.fai > ${GENOME_FILE}.chromsizes
fi

CHROM_SIZES=${GENOME_FILE}.chromsizes

density_file_basename=$(basename ${DENSITY_FILE})
density_file_basename=${density_file_basename[0]%%.*}

file_num=${#DENSITY_FILE[@]}

if [[ "${file_num}" == "1" ]]
then
    if [ ! -f "${genome_basename}.${WINDOWS}windows" ]; then
        bedtools makewindows -g ${CHROM_SIZES} -w ${WINDOWS} > ${genome_basename}.${WINDOWS}windows
    fi
    if [[ "${DENSITY_FILE[0]}" == "gc" ]]; then
        echo 'gc content'
        bedtools nuc -fi ${GENOME_FILE} -bed ${genome_basename}.${WINDOWS}windows > ${genome_basename}.${WINDOWS}windows.gc
        awk -v OFS='\t' '{ if (NR > 1) {print $1, $2, $3, $5}}' ${genome_basename}.${WINDOWS}windows.gc > ${genome_basename}.${WINDOWS}windows.gc.bg
        bedGraphToBigWig ${genome_basename}.${WINDOWS}windows.gc.bg ${CHROM_SIZES} ${density_file_basename}.bw
    else
        echo 'single file but not gc, treating input as bed'
        bedtools coverage -a ${genome_basename}.${WINDOWS}windows \
            -b ${DENSITY_FILE[0]} -g ${CHROM_SIZES} | awk -v OFS='\t' '{print $1, $2, $3, $4}' > ${density_file_basename}.bg

        bedGraphToBigWig ${density_file_basename}.bg ${CHROM_SIZES} ${density_file_basename}.bw
    fi
elif [[ "${file_num}" == "2" ]]
then
    echo 'two files, treating input as bam, where first is treatment, second is input'
    bamCompare -p ${THREADS} -b1 ${DENSITY_FILE[0]} -b2 ${DENSITY_FILE[1]} -o ${density_file_basename}.bw -of bigwig --scaleFactorsMethod readCount
else
    echo 'too many files inputted'
    exit 1
fi

if [[ "${SV}" == "TRA" ]]
then
computeMatrix reference-point -p ${THREADS} -S ${density_file_basename}.bw \
                            -R ${REGIONS_BED} \
                            --beforeRegionStartLength 10000 \
                            --referencePoint TSS \
                            --afterRegionStartLength 10000 \
                            -o ${OUTPUT_NAME}.mat.gz
else
computeMatrix scale-regions -p ${THREADS} -S ${density_file_basename}.bw \
                            -R ${REGIONS_BED} \
                            --beforeRegionStartLength 10000 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 10000 \
                            -o ${OUTPUT_NAME}.mat.gz
fi

plotProfile -m ${OUTPUT_NAME}.mat.gz \
            -out ${OUTPUT_NAME}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${OUTPUT_NAME}" \
            --outFileNameData ${OUTPUT_NAME}.tab