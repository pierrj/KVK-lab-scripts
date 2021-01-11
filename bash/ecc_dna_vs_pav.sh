while getopts g:s:a:t:l:f:e:m:n: option
do
case "${option}"
in
s) SAMPLE=${OPTARG};;
f) GFF_FILE=${OPTARG};;
e) ECCDNA_MAPFILE=${OPTARG};;
m) SAMPLE_MAPFILE=${OPTARG};;
n) ECC_NORMALIZATION=${OPTARG};; ## g for gene length multiplication or a for any overlap
p) PAV_FILE=${OPTARG};; ## file of genes (must match gene names in gff file) and number of genomes genes are found in
esac
done

awk '{if ($3 == "gene") print $0}' ${GFF_FILE} > ${basename_gff_file}.justgenes
awk '{if ($3 == "gene") print $0}' ${GFF_FILE} | awk -v OFS='\t' '{print substr($9,4, 10), $5-$4}' | awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | sort -k1,1 | awk '{print $2/1000}' > ${basename_gff_file}.gene_lengths

## look at confirmed spit reads per gene in all technical replicates
## normalize to limit bias against small genes which are more likely to be found in eccDNAs
if [ -f "${SAMPLE}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
fi
if [[ "${ECC_NORMALIZATION}" == "g" ]] ## normalize for gene length by multiplying by gene length 
then
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    bedtools intersect -f 1 -wa -c -a ${basename_gff_file}.justgenes -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $9, $10}' > ${ecc_basename}.splitreadspergene
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    paste ${ecc_basename}.splitreadspergene ${basename_gff_file}.gene_lengths | awk -v N=$num_srs '{print $1, ($2*$3)/N}' > ${ecc_basename}.normalized.splitreadspergene ## NORMALIZE TO DEAL WITH FAVORING OF SMALLER GENES TEST THIS LATER
    echo ${ecc_basename}.normalized.splitreadspergene >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE}
elif [[ "${ECC_NORMALIZATION}" == "a" ]] ## normalize for gene length by counting any overlap during bedtools intersect
then
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    bedtools intersect -wa -c -a ${basename_gff_file}.justgenes -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $9, $10}' > ${ecc_basename}.splitreadspergene
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    paste ${ecc_basename}.splitreadspergene ${basename_gff_file}.gene_lengths | awk -v N=$num_srs '{print $1, $2/N}' > ${ecc_basename}.normalized.splitreadspergene ## NORMALIZE TO DEAL WITH FAVORING OF SMALLER GENES TEST THIS LATER
    echo ${ecc_basename}.normalized.splitreadspergene >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE}
else
echo "invalid normalization type"
fi

# normalize and average across technical and biological replicates as written in previous scripts
if [ -f "${SAMPLE}.normalize_table_column" ]; then
    rm ${SAMPLE}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${SAMPLE}.normalize_table_column ; done
paste ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn ${SAMPLE}.normalize_table_column ${SAMPLE_MAPFILE} > ${SAMPLE}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${SAMPLE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.normalized.splitreadspergene

## generate output file
# gene splitreads versus PAV (PAV FILE MUST BE IN EXACT SAME ORDER AS GFF FILE)
awk '{print $2}' ${SAMPLE}.normalized.splitreadspergene > ${SAMPLE}.normalized.splitreadspergene.countcolumn
paste ${PAV_FILE} ${SAMPLE}.normalized.splitreadspergene.countcolumn > ${SAMPLE}.PAVvsSRs

Rscript --vanilla /global/home/users/pierrj/git/R/basic_scatterplot.R ${SAMPLE}.PAVvsSRs 2 3 PAV SRs ${SAMPLE}.PAVvsSRs ${SAMPLE}.PAVvsSRs 4 4