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
m) SAMPLE_MAPFILE=${OPTARG};;
esac
done

# make 100kb bins and count repeats per bin
samtools faidx ${GENOME_FASTA}
cut -f1,2 ${GENOME_FASTA}.fai > ${genome_fasta_basename}.sizes
bedtools makewindows -g ${genome_fasta_basename}.sizes -w 100000 | awk '$3-$2==100000' > ${genome_fasta_basename}.100kbins # no bins smaller than 100kb
bedtools intersect -a ${genome_fasta_basename}.100kbins -b ${REPEAT_FILE} -c > ${SAMPLE}.repeatsper100kb

if [ -f "${SAMPLE}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    bedtools intersect -a ${genome_fasta_basename}.100kbins -b ${ECCDNA_FILE} -c > ${ecc_basename}.eccsper100kb
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    awk -v N=$num_srs '{print $1, $2, $4/N}' ${ecc_basename}.eccsper100kb > ${ecc_basename}.eccsper100kb.normalized
    echo ${ecc_basename}.eccsper100kb.normalized >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE}

# normalize and average across technical and biological replicates as written in previous scripts
if [ -f "${SAMPLE}.normalize_table_column" ]; then
    rm ${SAMPLE}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${SAMPLE}.normalize_table_column ; done
paste ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn ${SAMPLE}.normalize_table_column ${SAMPLE_MAPFILE} > ${SAMPLE}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${SAMPLE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.normalized.splitreadsper100kb

# look at scaffold averages instead of 100kb bins
awk -v OFS='\t' '{seen[$1]+=$4; count[$1]++} END{for (x in seen)print x, seen[x]/count[x]}' ${SAMPLE}.repeatsper100kb | sort -k1,1 > ${SAMPLE}.repeatsper100kb.scaffoldaverage
awk -v OFS='\t' '{seen[$1]+=$3; count[$1]++} END{for (x in seen)print x, seen[x]/count[x]}' ${SAMPLE}.normalized.splitreadsper100kb | sort -k1,1 > ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage

# number of eccDNA forming regions vs repeats per 100kb bins 
awk '{print $3}' ${SAMPLE}.normalized.splitreadsper100kb > ${SAMPLE}.normalized.splitreadsper100kb.countcolumn
paste ${SAMPLE}.repeatsper100kb ${SAMPLE}.normalized.splitreadsper100kb.countcolumn > ${SAMPLE}.SRsvsrepeatsper100kb
# number of repeats per 100kb bins versus eccDNA forming regions (averages per scaffold)
awk '{print $2}' ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage > ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage.countcolumn
paste ${SAMPLE}.repeatsper100kb.scaffoldaverage ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage.countcolumn > ${SAMPLE}.SRsvsrepeatsper100kbperscaffold

Rscript --vanilla /global/home/users/pierrj/git/R/basic_scatterplot.R ${SAMPLE}.SRsvsrepeatsper100kb 4 5 repeatsper100kb SRs ${SAMPLE}.SRsvsrepeatsper100kb ${SAMPLE}.SRsvsrepeatsper100kb 4 4

Rscript --vanilla /global/home/users/pierrj/git/R/basic_scatterplot.R ${SAMPLE}.SRsvsrepeatsper100kbperscaffold 2 3 RepeatsPerScaffold SRs ${SAMPLE}.SRsvsrepeatsper100kbperscaffold ${SAMPLE}.SRsvsrepeatsper100kbperscaffold 4 4 