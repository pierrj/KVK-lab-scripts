#!/bin/bash
while getopts g:s:a:t:l:f:e:m: option
do
case "${option}"
in
g) GENOME_FASTA=${OPTARG};;
s) SAMPLE=${OPTARG};;
a) SRA_LIST=${OPTARG};;
t) THREADS=${OPTARG};;
l) LIBTYPE=${OPTARG};; ## 1 for SE or 2 for PE
f) GFF_FILE=${OPTARG};;
e) ECCDNA_MAPFILE=${OPTARG};;
m) SAMPLE_MAPFILE=${OPTARG};;
esac
done

genome_fasta_basename=$(basename ${GENOME_FASTA})

if [ -d "${genome_fasta_basename}_starindex" ]; then
    rm -r ${genome_fasta_basename}_starindex
fi

mkdir ${genome_fasta_basename}_starindex

## index genome for STAR
STAR --runThreadN ${THREADS} --runMode genomeGenerate --genomeDir ${genome_fasta_basename}_starindex \
    --genomeFastaFiles ${GENOME_FASTA} \
    --sjdbGTFfile ${GFF_FILE} \
    --sjdbOverhang 100 \
    --genomeSAindexNbases 11 \
    --sjdbGTFtagExonParentTranscript ID \
    --sjdbGTFtagExonParentGene Parent

## get exon and gene length per gene (ONLY WORKS IF YOU HAVE ONE TRANSCRIPT PER GENE AND SPECIFIC FORMAT FROM FUNGAP)
basename_gff_file=$(basename ${GFF_FILE})
grep 'exon' ${GFF_FILE} | awk -v OFS='\t' '{print substr($9,4, 10), $5-$4}' | awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | sort -k1,1 | awk '{print $2/1000}' > ${basename_gff_file}.exon_lengths
grep 'gene' ${GFF_FILE} > ${basename_gff_file}.justgenes
grep 'gene' ${GFF_FILE} | awk -v OFS='\t' '{print substr($9,4, 10), $5-$4}' | awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | sort -k1,1 | awk '{print $2/1000}' > ${basename_gff_file}.gene_lengths

## download and map all reads in passed list of SRA accessions
## generate RPKM for each gene
while read SRA; do
    if [[ "${LIBTYPE}" -eq 1 ]]
    then
        /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
        /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
        STAR --runThreadN ${THREADS} \
            --genomeDir ${genome_fasta_basename}_starindex \
            --readFilesIn ${SRA}.sra.fastq \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${SRA}. \
            --quantMode GeneCounts
    elif [[ "${LIBTYPE}" -eq 2 ]]
    then
        /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SRA} -O .
        /global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SRA}.sra
        STAR --runThreadN ${THREADS} \
            --genomeDir ${genome_fasta_basename}_starindex \
            --readFilesIn ${SRA}.sra_1.fastq ${SRA}.sra_2.fastq \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${SRA}. \
            --quantMode GeneCounts
    else
    echo "invalid libtype"
    fi
    awk 'NR>4 {print $2}' ${SRA}.ReadsPerGene.out.tab > ${SRA}.ReadsPerGene.out.genecolumn.tab
    num_reads=$(awk '{SUM+=$1}END{print SUM/1000000}' ${SRA}.ReadsPerGene.out.genecolumn.tab)
    paste ${SRA}.ReadsPerGene.out.genecolumn.tab ${basename_gff_file}.exon_lengths | awk -v N=$num_reads '{print $1/($2*N)}' > ${SAMPLE}.RPKM.${SRA}.ReadsPerGene.out.genecolumn.tab
done < ${SRA_LIST}

## average RPKMs across input SRA accessions
first_SRA=$(head -1 ${SRA_LIST})
awk 'NR>4 {print $1}' ${first_SRA}.ReadsPerGene.out.tab > ${SAMPLE}.genecount_firstcolumn
paste $(find . -maxdepth 1 -name "${SAMPLE}.RPKM.*.ReadsPerGene.out.genecolumn.tab" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") > ${SAMPLE}.genecount_table
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' ${SAMPLE}.genecount_table > ${SAMPLE}.genecount_table_average
paste ${SAMPLE}.genecount_firstcolumn ${SAMPLE}.genecount_table_average_lengthnormalized > ${SAMPLE}.genecount_table_final

## look at confirmed spit reads per gene in all technical replicates
## normalize to limit bias against small genes which are more likely to be found in eccDNAs
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    bedtools intersect -f 1 -wa -c -a ${basename_gff_file}.justgenes -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $9, $10}' > ${ecc_basename}.splitreadspergene #### CHECK THE COLUMNS HERE, should be gene name and count per gene
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    paste ${ecc_basename}.splitreadspergene ${basename_gff_file}.gene_lengths | awk -v N=$num_srs '{print $1, ($2*$3)/N)}' > ${ecc_basename}.normalized.splitreadspergene ## NORMALIZE TO DEAL WITH FAVORING OF SMALLER GENES DOUBLE CHECK THIS
done < ${ECCDNA_MAPFILE}
# normalize and average across technical and biological replicates as written in previous scripts
if [ -f "${SAMPLE}.normalize_table" ]; then
    rm ${SAMPLE}.normalize_table
fi
sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk {print $1+1})
for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${SAMPLE}.normalize_table ; done
ECCDNA_FILE=$(head -1 ${ECCDNA_MAPFILE})
ecc_basename=$(basename ${ECCDNA_FILE})
/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -t ${ecc_basename}.normalized.splitreadspergene -m ${SAMPLE_MAPFILE} -n ${SAMPLE}.normalize_table -y t
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.normalized.splitreadspergene

# make 100kb bins and count genes per bin
samtools faidx ${GENOME_FASTA}
cut -f1,2 ${GENOME_FASTA}.fai > ${genome_fasta_basename}.sizes
bedtools makewindows -g ${genome_fasta_basename}.sizes -w 100000 | awk '$3-$2==100000' > ${genome_fasta_basename}.100kbins # no bins smaller than 100kb
bedtools intersect -a ${genome_fasta_basename}.100kbins -b ${basename_gff_file}.justgenes -c > ${SAMPLE}.genesperk100kb

while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    bedtools intersect -a ${genome_fasta_basename}.100kbbins -b ${ECCDNA_FILE} -c > ${ecc_basename}.eccsper100kb
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    awk -v N=$num_srs '{print $1, $2, $3/N}' > ${ecc_basename}.eccsper100kb.normalized ## DOUBLE CHECK COLUMN COUNT HERE
done < ${ECCDNA_MAPFILE}
ECCDNA_FILE=$(head -1 ${ECCDNA_MAPFILE})
ecc_basename=$(basename ${ECCDNA_FILE})
/global/home/users/pierrj/git/bash/create_mapfile_for_normalize_and_average.sh -t ${ecc_basename}.eccsper100kb.normalized -m ${SAMPLE_MAPFILE} -n ${SAMPLE}.normalize_table -y t
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.normalized.splitreadsper100kb

# look at scaffold averages instead of 100kb bins
awk '{seen[$1]+=$4; count[$1]++} END{for (x in seen)print x, seen[x]/count[x]}' ${SAMPLE}.genesperk100kb | sort -k1,1 > ${SAMPLE}.genesper100kb.scaffoldaverage
awk '{seen[$1]+=$4; count[$1]++} END{for (x in seen)print x, seen[x]/count[x]}' ${SAMPLE}.normalized.splitreadsper100kb | sort -k1,1 > ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage

## generate output files
# gene splitreads versus RPKM
awk '{print $2}' ${SAMPLE}.normalized.splitreadspergene > ${SAMPLE}.normalized.splitreadspergene.countcolumn
paste ${SAMPLE}.genecount_table_final ${SAMPLE}.normalized.splitreadspergene.countcolumn > ${SAMPLE}.RPKMvsSRs
# number of eccDNA forming regions vs genes per 100kb bins 
awk '{print $3}' ${SAMPLE}.normalized.splitreadsper100kb > ${SAMPLE}.normalized.splitreadsper100kb.countcolumn
paste ${SAMPLE}.genesperk100kb ${SAMPLE}.normalized.splitreadsper100kb.countcolumn > ${SAMPLE}.SRsvsgenesper100kb
# number of genes per 100kb bins versus eccDNA forming regions (averages per scaffold)
awk '{print $4}' ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage > ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage.countcolumn
paste ${SAMPLE}.genesper100kb.scaffoldaverage ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage.countcolumn > ${SAMPLE}.SRsvsgenesperk100kbperscaffold


## DOUBLE CHECK THAT GENE EXON LENGTHS AND READS PER GENE LINE UP PERFECTLY