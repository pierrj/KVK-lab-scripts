#!/bin/bash
while getopts o:g:m: option
do
case "${option}"
in
o) OUTPUT=${OPTARG};;
g) GENOME=${OPTARG};;
m) MAPFILE=${OPTARG};;
esac
done

if [ -f "${OUTPUT}" ]; then
    rm ${OUTPUT}
fi

while read sample;
do
    cd ${sample}
    ecc_count=$(wc -l ecccaller_output.${sample}.renamed.details.tsv | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${GENOME} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${OUTPUT}
    cd ..
done < ${MAPFILE}