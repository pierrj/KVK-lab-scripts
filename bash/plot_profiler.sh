#!/bin/bash
while getopts b:d:g:w:t:o: option
do
case "${option}"
in
b) REGIONS_BED=${OPTARG};;
d) DENSITY_BED=${OPTARG};;
g) GENOME_FILE=${OPTARG};;
w) WINDOWS=${OPTARG};;
t) THREADS=${OPTARG};;
o) OUTPUT_NAME=${OPTARG};;
esac
done


genome_basename=$(basename ${GENOME_FILE})

if [ ! -f "${GENOME_FILE}.chromsizes" ]; then
    samtools faidx ${GENOME_FILE}
    cut -f1,2 ${GENOME_FILE}.fai > ${GENOME_FILE}.chromsizes
fi

CHROM_SIZES=${GENOME_FILE}.chromsizes

if [ ! -f "${genome_basename}.${WINDOWS}windows" ]; then
    bedtools makewindows -g ${CHROM_SIZES} -w ${WINDOWS} > ${genome_basename}.${WINDOWS}windows
fi


## make windows first

density_bed_basename=$(basename ${DENSITY_BED} )


bedtools coverage -a ${genome_basename}.${WINDOWS}windows \
    -b ${DENSITY_BED} -g ${CHROM_SIZES} | awk -v OFS='\t' '{print $1, $2, $3, $4}' > ${density_bed_basename}.bg

bedGraphToBigWig ${density_bed_basename}.bg ${CHROM_SIZES} ${density_bed_basename}.bw


## SKIP ZEROS OR NO?
computeMatrix scale-regions -p ${THREADS} -S ${density_bed_basename}.bw \
                            -R ${REGIONS_BED} \
                            --beforeRegionStartLength 10000 \
                            --regionBodyLength 1000 \
                            --afterRegionStartLength 10000 \
                            -o ${OUTPUT_NAME}.mat.gz

plotProfile -m ${OUTPUT_NAME}.mat.gz \
            -out ${OUTPUT_NAME}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${OUTPUT_NAME}" \
            --outFileNameData ${OUTPUT_NAME}.tab
