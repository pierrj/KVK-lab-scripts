#!/bin/bash
while getopts g:s:a:t:f:e:m: option
do
case "${option}"
in
g) GENOME_FASTA=${OPTARG};;
s) SAMPLE=${OPTARG};;
a) SRA_LIST=${OPTARG};;
t) THREADS=${OPTARG};;
f) REPEAT_FILE=${OPTARG};;
e) ECCDNA_MAPFILE=${OPTARG};;
d) ECCDNA_MAPFILE_DETAILS=${OPTARG};;
m) SAMPLE_MAPFILE=${OPTARG};;
esac
done

basename_gff_file=$(basename ${GFF_FILE})
awk '{if ($3 == "gene") print $0}' ${GFF_FILE} > ${basename_gff_file}.justgenes


## get ecc density
if [ -f "${SAMPLE}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE BAMFILE; do
    basename_eccdna_file=$(basename ${ECCDNA_FILE})
    if [ -f "${basename_eccdna_file}.ecc_density" ]; then
        rm ${basename_eccdna_file}.ecc_density
    fi
    bedtools intersect -u -f 0.5 -a ${ECCDNA_FILE} -b ${basename_gff_file}.justgenes > genic.${basename_eccdna_file}
    bedtools intersect -v -f 0.5 -a ${ECCDNA_FILE} -b ${basename_gff_file}.justgenes > noncoding.${basename_eccdna_file}
    ecc_count_genic=$(awk ' $6 != "lowq" ' genic.${basename_eccdna_file} | wc -l | awk '{print $1}')
    ecc_count_noncoding=$(awk ' $6 != "lowq" ' noncoding.${basename_eccdna_file} | wc -l | awk '{print $1}')
    ecc_count_all=$(awk ' $6 != "lowq" ' ${ECCDNA_FILE} | wc -l | awk '{print $1}')
    read_count=$(samtools view -c ${BAMFILE} | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${GENOME_FASTA} ) ## megabasepairs
    total_genic=$(awk -v ECC="$ecc_count_genic" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    total_noncoding=$(awk -v ECC="$ecc_count_noncoding" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    total_all=$(awk -v ECC="$ecc_count_all" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepaircd
    echo -e genic'\t'$total_genic >> ${basename_eccdna_file}.ecc_density
    echo -e noncoding'\t'$total_noncoding >> ${basename_eccdna_file}.ecc_density
    echo -e all'\t'$total_all >> ${basename_eccdna_file}.ecc_density
    echo ${basename_eccdna_file}.ecc_density >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE_DETAILS}

## average across tech reps and bioreps
if [ -f "${SAMPLE}.normalize_table_column" ]; then
    rm ${SAMPLE}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${SAMPLE}.normalize_table_column ; done
paste ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn ${SAMPLE}.normalize_table_column ${SAMPLE_MAPFILE} > ${SAMPLE}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${SAMPLE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.ecc_density

## get average length of eccDNAs (weighted by frequency)
if [ -f "${SAMPLE}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE; do
    basename_eccdna_file=$(basename ${ECCDNA_FILE})
    average_length=$(awk '{print $3-$2}' ${ECCDNA_FILE} | awk '{sum+=$1} END {print sum/NR}')
    echo -e average_length_weighted'\t'${basename_eccdna_file}.average_length_weighted
    echo ${basename_eccdna_file}.average_length_weighted >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE_SRS}
if [ -f "${SAMPLE}.normalize_table_column" ]; then
    rm ${SAMPLE}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${SAMPLE}.normalize_table_column ; done
paste ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn ${SAMPLE}.normalize_table_column ${SAMPLE_MAPFILE} > ${SAMPLE}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${SAMPLE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.ecc_density


## get average length of eccDNAs (unweighted)
if [ -f "${SAMPLE}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE BAMFILE; do
    basename_eccdna_file=$(basename ${ECCDNA_FILE})
    average_length=$(awk '{print $3-$2}' ${ECCDNA_FILE} | awk '{sum+=$1} END {print sum/NR}')
    echo -e average_length_weighted'\t'${basename_eccdna_file}.average_length_weighted
    echo ${basename_eccdna_file}.average_length_weighted >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE_DETAILS}
if [ -f "${SAMPLE}.normalize_table_column" ]; then
    rm ${SAMPLE}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${SAMPLE}.normalize_table_column ; done
paste ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn ${SAMPLE}.normalize_table_column ${SAMPLE_MAPFILE} > ${SAMPLE}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${SAMPLE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.ecc_density

## get length distribution histogram

## determine max_max_length
if [ -f "${SAMPLE}.ecc_maxlength" ]; then
    rm ${SAMPLE}.ecc_maxlength
fi
while read ECCDNA_FILE; do
    basename_eccdna_file=$(basename ${ECCDNA_FILE})
    max_length=$(awk '{print $3-$2}' ${ECCDNA_FILE} | sort -k1,1nr | head -1)
    echo max_length >> ${SAMPLE}.ecc_maxlength
done < ${ECCDNA_MAPFILE_SRS}
bin_size=100
max_max_length=$(sort -k1,1nr ${SAMPLE}.ecc_maxlength | head -1 | awk -v b="$bin_size" '{print $0/b}' | awk '{print ($0-int($0)>0)?int($0)+1:int($0)}' | awk -v b="$bin_size" '{print $0*b}')

if [ -f "${SAMPLE}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE; do
    basename_eccdna_file=$(basename ${ECCDNA_FILE})
    awk '{print $3-$2}' ${ECCDNA_FILE} > ${basename_eccdna_file}.ecc_lengths
    python /global/home/users/pierrj/git/python/get_density.py ${basename_eccdna_file}.ecc_lengths $bin_size $max_max_length ${basename_eccdna_file}.ecc_lengths_density ## PATH NEEDS TO BE FIXED
    echo ${basename_eccdna_file}.ecc_lengths_density >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE_SRS}
if [ -f "${SAMPLE}.normalize_table_column" ]; then
    rm ${SAMPLE}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${SAMPLE}.normalize_table_column ; done
paste ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn ${SAMPLE}.normalize_table_column ${SAMPLE_MAPFILE} > ${SAMPLE}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${SAMPLE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.ecc_density

## get flanking repeats stuff