
## FOR CT

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/gene_analysis_wrapper/CT/raw_files/test

treatment=CT
MAPFILE=/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/mapfile_${treatment}
GENE_BEDFILE=/global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.just_genes.post_repeatmasking.T0name.bed
OUTPUT_NAME=CT
PANNZER_INPUT=/global/scratch/users/pierrj/references/pannzer_GUY11_PacBio_merge_KROJ_post_repeatmasking.txt
SCORE_CUTOFF=0.6

rm *

while read SAMPLE; do cp /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/${SAMPLE}/ecccaller_output.${SAMPLE}.renamed.details.tsv . ; done < ${MAPFILE}

while read SAMPLE
do
    #awk '$3-$2>1000' ecccaller_output.${SAMPLE}.renamed.details.tsv > ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv
    cat ecccaller_output.${SAMPLE}.renamed.details.tsv > ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv
    line_count=$(wc -l ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv | awk '{print int($1/100)}') 
    sr_count_cutoff=$(sort -k4,4nr ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv | head -${line_count} | tail -1 | awk '{print $4}')
    awk -v CUTOFF=${sr_count_cutoff} '$4>=CUTOFF' ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv > ecccaller_output.${SAMPLE}.common.renamed.details.tsv
    #cat ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv > ecccaller_output.${SAMPLE}.common.renamed.details.tsv
done < ${MAPFILE}

## TO GENERATE BIOREP MAPFILE IF STANDARD FILES ARE SET

awk '{print substr($1, 0,4)}' ${MAPFILE} | sort | uniq > tmp_biorepmapfile

while read bio_rep; 
do
    cat $(find . -maxdepth 1 -name "ecccaller_output.${bio_rep}*.common.renamed.details.tsv" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") > ecccaller_output.${bio_rep}.common.renamed.details.tsv
done < tmp_biorepmapfile

while read bio_rep; 
do
    bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b ecccaller_output.${bio_rep}.common.renamed.details.tsv | awk '{if ($5!=0) {print $4}}' > ${bio_rep}.common.genes
done < tmp_biorepmapfile

cat $(find . -maxdepth 1 -name "*.common.genes" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | sort | uniq -c | awk '$1==3 {print $2}' > ${OUTPUT_NAME}.common.genes
#cat $(find . -maxdepth 1 -name "*.common.genes" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | sort | uniq -c | awk '{if ($1==3 || $1==2) {print $2}}' > ${OUTPUT_NAME}.common.genes


## FOR PQ
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/gene_analysis_wrapper/PQ/raw_files/test


treatment=PQ
MAPFILE=/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/mapfile_${treatment}
GENE_BEDFILE=/global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.just_genes.post_repeatmasking.T0name.bed
OUTPUT_NAME=PQ
PANNZER_INPUT=/global/scratch/users/pierrj/references/pannzer_GUY11_PacBio_merge_KROJ_post_repeatmasking.txt
SCORE_CUTOFF=0.6

rm *

while read SAMPLE; do cp /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/${SAMPLE}/ecccaller_output.${SAMPLE}.renamed.details.tsv . ; done < ${MAPFILE}

while read SAMPLE
do
    #awk '$3-$2>1000' ecccaller_output.${SAMPLE}.renamed.details.tsv > ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv
    cat ecccaller_output.${SAMPLE}.renamed.details.tsv > ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv
    line_count=$(wc -l ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv | awk '{print int($1/100)}') 
    sr_count_cutoff=$(sort -k4,4nr ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv | head -${line_count} | tail -1 | awk '{print $4}')
    awk -v CUTOFF=${sr_count_cutoff} '$4>=CUTOFF' ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv > ecccaller_output.${SAMPLE}.common.renamed.details.tsv
    #cat ecccaller_output.${SAMPLE}.hconf.renamed.details.tsv > ecccaller_output.${SAMPLE}.common.renamed.details.tsv
done < ${MAPFILE}

## TO GENERATE BIOREP MAPFILE IF STANDARD FILES ARE SET

awk '{print substr($1, 0,4)}' ${MAPFILE} | sort | uniq > tmp_biorepmapfile

while read bio_rep; 
do
    cat $(find . -maxdepth 1 -name "ecccaller_output.${bio_rep}*.common.renamed.details.tsv" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") > ecccaller_output.${bio_rep}.common.renamed.details.tsv
done < tmp_biorepmapfile

while read bio_rep; 
do
    bedtools intersect -f 1 -wa -c -a ${GENE_BEDFILE} -b ecccaller_output.${bio_rep}.common.renamed.details.tsv | awk '{if ($5!=0) {print $4}}' > ${bio_rep}.common.genes
done < tmp_biorepmapfile

cat $(find . -maxdepth 1 -name "*.common.genes" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | sort | uniq -c | awk '$1==3 {print $2}' > ${OUTPUT_NAME}.common.genes

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ctvspqgenecomparison

ECCDNA_MAPFILE=CT.eccdna_mapfile
SAMPLE=CT
SAMPLE_MAPFILE=CT.sample_mapfile
GFF_FILE=sorted.ctvspqhotgene.gff3

if [ -f "${SAMPLE}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    bedtools intersect -f 1 -wa -c -a ${GFF_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $9, $10}' > ${ecc_basename}.splitreadspergene
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    awk -v N=$num_srs '{print $1, $2/N}' ${ecc_basename}.splitreadspergene > ${ecc_basename}.normalized.splitreadspergene
    echo ${ecc_basename}.normalized.splitreadspergene >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE}

# normalize and average across technical and biological replicates as written in previous scripts
if [ -f "${SAMPLE}.normalize_table_column" ]; then
    rm ${SAMPLE}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${SAMPLE}.normalize_table_column ; done
paste ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn ${SAMPLE}.normalize_table_column ${SAMPLE_MAPFILE} > ${SAMPLE}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${SAMPLE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.normalized.splitreadspergene

ECCDNA_MAPFILE=PQ.eccdna_mapfile
SAMPLE=PQ
SAMPLE_MAPFILE=PQ.sample_mapfile
GFF_FILE=sorted.ctvspqhotgene.gff3

if [ -f "${SAMPLE}.mapfile_for_normalize_and_average_filecolumn" ]; then
    rm ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
fi
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    bedtools intersect -f 1 -wa -c -a ${GFF_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $9, $10}' > ${ecc_basename}.splitreadspergene
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    awk -v N=$num_srs '{print $1, $2/N}' ${ecc_basename}.splitreadspergene > ${ecc_basename}.normalized.splitreadspergene
    echo ${ecc_basename}.normalized.splitreadspergene >> ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn
done < ${ECCDNA_MAPFILE}

# normalize and average across technical and biological replicates as written in previous scripts
if [ -f "${SAMPLE}.normalize_table_column" ]; then
    rm ${SAMPLE}.normalize_table_column
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${SAMPLE}.normalize_table_column ; done
paste ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn ${SAMPLE}.normalize_table_column ${SAMPLE_MAPFILE} > ${SAMPLE}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${SAMPLE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.normalized.splitreadspergene

awk -v OFS='\t' '{print $2, $3}' PQ.normalized.splitreadspergene > PQ.normalized.splitreadspergene.justaverageandsd
paste CT.normalized.splitreadspergene PQ.normalized.splitreadspergene.justaverageandsd > srspergene_ctvspq

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
genome_fasta_basename=$(basename ${GENOME_FASTA})
GFF_FILE=sorted.ctvspqhotgene.gff3

samtools faidx ${GENOME_FASTA}
cut -f1,2 ${GENOME_FASTA}.fai > ${genome_fasta_basename}.sizes
bedtools makewindows -g ${genome_fasta_basename}.sizes -w 100000 | awk '$3-$2==100000' > ${genome_fasta_basename}.100kbins

bedtools intersect -a ${genome_fasta_basename}.100kbins -b ${GFF_FILE} -c > hotgenesper100kb