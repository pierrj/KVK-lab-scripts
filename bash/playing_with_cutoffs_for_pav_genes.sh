while read gene; do
grep ${gene} G3_1A.PAVvsSRs_aoption >> G3_1A.PAVvsSRs_aoption.effectorsonly
done < guy11_effector_protein_names

if [ -f "${SAMPLE}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    bedtools intersect -f 1 -wa -c -a ${basename_gff_file}.justgenes -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $9, $10}' > ${ecc_basename}.splitreadspergene
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    paste ${ecc_basename}.splitreadspergene ${basename_gff_file}.gene_lengths | awk -v N=$num_srs '{print $1, $2/N}' > ${ecc_basename}.normalized.splitreadspergene ## NORMALIZE TO DEAL WITH FAVORING OF SMALLER GENES TEST THIS LATER
    echo ${ecc_basename}.normalized.splitreadspergene >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE}


if [ -f "${SAMPLE}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    awk '$3-$2 > 1000' ${ECCDNA_FILE} > lengthfiltered.${ecc_basename}
    bedtools intersect -wa -c -a ${basename_gff_file}.justgenes -b lengthfiltered.${ecc_basename} | awk -v OFS='\t' '{print $9, $10}' > ${ecc_basename}.splitreadspergene
    num_srs=$(wc -l lengthfiltered.${ecc_basename} | awk '{print $1/100000}')
    paste ${ecc_basename}.splitreadspergene ${basename_gff_file}.gene_lengths | awk -v N=$num_srs '{print $1, $2/N}' > ${ecc_basename}.normalized.splitreadspergene ## NORMALIZE TO DEAL WITH FAVORING OF SMALLER GENES TEST THIS LATER
    echo ${ecc_basename}.normalized.splitreadspergene >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE}