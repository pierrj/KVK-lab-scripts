#!/bin/bash
while getopts b:d:g:w:t:o:s:f: option
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
f) FLANKING_DIST=${OPTARG};;
esac
done

genome_basename=$(basename ${GENOME_FILE})
genome_basename=${genome_basename%%.*}

samtools faidx ${GENOME_FILE}
cut -f1,2 ${GENOME_FILE}.fai > ${GENOME_FILE}.chromsizes

CHROM_SIZES=${GENOME_FILE}.chromsizes

density_file_basename=$(basename ${DENSITY_FILE})
density_file_basename=${density_file_basename[0]%%.*}

file_num=${#DENSITY_FILE[@]}

bedtools makewindows -g ${CHROM_SIZES} -w ${WINDOWS} > ${genome_basename}.${WINDOWS}windows

if [[ "${file_num}" == "1" ]]
then
    if [[ "${DENSITY_FILE[0]}" == "gc" ]]; then
        echo 'gc content'
        bedtools nuc -fi ${GENOME_FILE} -bed ${genome_basename}.${WINDOWS}windows > ${genome_basename}.${WINDOWS}windows.gc
        ## THIS IS ACTUALLY GC CONTENT ##
        awk -v OFS='\t' '{ if (NR > 1) {print $1, $2, $3, $5}}' ${genome_basename}.${WINDOWS}windows.gc > ${genome_basename}.${WINDOWS}windows.gc.bg
        bedGraphToBigWig ${genome_basename}.${WINDOWS}windows.gc.bg ${CHROM_SIZES} ${density_file_basename}.bw
    elif [[ "${DENSITY_FILE[0]}" == *"bismark"* ]]; then
        echo 'methylation'
        bedGraphToBigWig ${DENSITY_FILE[0]} ${CHROM_SIZES} ${density_file_basename}.bw
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
    echo 'more than two files, averaging all of them together'
    for density_file in "${DENSITY_FILE[@]}"; do
        basename_density_file=$(basename $density_file)
        read_count=$(samtools view -c -F 4 -F 2048 $density_file | awk '{print $1/1000000}')
        bedtools coverage -counts -sorted -a ${genome_basename}.${WINDOWS}windows \
            -b $density_file -g ${CHROM_SIZES} | awk -v r=$read_count -v OFS='\t' '{print $(NF)/r}' > ${basename_density_file}.${OUTPUT_NAME}.bg
    done
    paste *.${OUTPUT_NAME}.bg | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > averaged.${OUTPUT_NAME}.bg
    paste ${genome_basename}.${WINDOWS}windows averaged.${OUTPUT_NAME}.bg > ${density_file_basename}.bg
    bedGraphToBigWig ${density_file_basename}.bg ${CHROM_SIZES} ${density_file_basename}.bw
fi

if [[ "${SV}" == "TRA" ]]
then
computeMatrix reference-point -p ${THREADS} -S ${density_file_basename}.bw \
                            -R ${REGIONS_BED} \
                            --beforeRegionStartLength $FLANKING_DIST \
                            --referencePoint TSS \
                            --afterRegionStartLength $FLANKING_DIST \
                            -o ${OUTPUT_NAME}.mat.gz
else
computeMatrix scale-regions -p ${THREADS} -S ${density_file_basename}.bw \
                            -R ${REGIONS_BED} \
                            --beforeRegionStartLength $FLANKING_DIST \
                            --regionBodyLength $(( FLANKING_DIST / 2)) \
                            --afterRegionStartLength $FLANKING_DIST \
                            -o ${OUTPUT_NAME}.mat.gz
fi

plotProfile -m ${OUTPUT_NAME}.mat.gz \
            -out ${OUTPUT_NAME}.pdf \
            --numPlotsPerRow 1 \
            --plotTitle "${OUTPUT_NAME}" \
            --outFileNameData ${OUTPUT_NAME}.tab