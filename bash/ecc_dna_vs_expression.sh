
#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.#!/bin/bash
while getopts g:s:a:t:l:f:e:m:n: option
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
n) ECC_NORMALIZATION=${OPTARG};; ## g for gene length multiplication or a for any overlap or n for no normalization
esac
done

## USAGE ##
# this script maps rnaseq data and compares it to junction split reads associated with eccdnas
# -g genome fasta file
# -s output short name
# -t threads for STAR
# -l library type, either single end or paired end
# -f gff file of gene annotations
# -e list of eccdna files
# -m list of eccdna sequencing sample names
# -n type of normalization for number of junction split reads for genes, usually no normalization makes sense

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
awk '{if ($3 == "gene") print $0}' ${GFF_FILE} > ${basename_gff_file}.justgenes
awk '{if ($3 == "gene") print $0}' ${GFF_FILE} | awk -v OFS='\t' '{print substr($9,4, 10), $5-$4}' | awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' | sort -k1,1 | awk '{print $2/1000}' > ${basename_gff_file}.gene_lengths


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
paste ${SAMPLE}.genecount_firstcolumn ${SAMPLE}.genecount_table_average > ${SAMPLE}.genecount_table_final

## look at junction split reads spit reads per gene in all technical replicates
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
elif [[ "${ECC_NORMALIZATION}" == "n" ]] ## dont normalize, only count full overlaps
then
while read ECCDNA_FILE; do
    ecc_basename=$(basename ${ECCDNA_FILE})
    bedtools intersect -f 1 -wa -c -a ${basename_gff_file}.justgenes -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $9, $10}' > ${ecc_basename}.splitreadspergene
    num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
    paste ${ecc_basename}.splitreadspergene ${basename_gff_file}.gene_lengths | awk -v N=$num_srs '{print $1, $2/N}' > ${ecc_basename}.normalized.splitreadspergene
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

# make 100kb bins and count genes per bin
samtools faidx ${GENOME_FASTA}
cut -f1,2 ${GENOME_FASTA}.fai > ${genome_fasta_basename}.sizes
bedtools makewindows -g ${genome_fasta_basename}.sizes -w 100000 | awk '$3-$2==100000' > ${genome_fasta_basename}.100kbins # no bins smaller than 100kb
bedtools intersect -a ${genome_fasta_basename}.100kbins -b ${basename_gff_file}.justgenes -c > ${SAMPLE}.genesperk100kb

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
paste ${SAMPLE}.mapfile_for_normalize_and_average_filecolumn ${SAMPLE}.normalize_table_column ${SAMPLE_MAPFILE} > ${SAMPLE}.mapfile_for_normalize_and_average
/global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${SAMPLE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 3 -n n
mv ${SAMPLE}.normalized_binned ${SAMPLE}.normalized.splitreadsper100kb

# look at scaffold averages instead of 100kb bins
awk -v OFS='\t' '{seen[$1]+=$4; count[$1]++} END{for (x in seen)print x, seen[x]/count[x]}' ${SAMPLE}.genesperk100kb | sort -k1,1 > ${SAMPLE}.genesper100kb.scaffoldaverage
awk -v OFS='\t' '{seen[$1]+=$3; count[$1]++} END{for (x in seen)print x, seen[x]/count[x]}' ${SAMPLE}.normalized.splitreadsper100kb | sort -k1,1 > ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage

## generate output files
# gene splitreads versus RPKM
awk '{print $2}' ${SAMPLE}.normalized.splitreadspergene > ${SAMPLE}.normalized.splitreadspergene.countcolumn
paste ${SAMPLE}.genecount_table_final ${SAMPLE}.normalized.splitreadspergene.countcolumn > ${SAMPLE}.RPKMvsSRs
# number of eccDNA forming regions vs genes per 100kb bins 
awk '{print $3}' ${SAMPLE}.normalized.splitreadsper100kb > ${SAMPLE}.normalized.splitreadsper100kb.countcolumn
paste ${SAMPLE}.genesperk100kb ${SAMPLE}.normalized.splitreadsper100kb.countcolumn > ${SAMPLE}.SRsvsgenesper100kb
# number of genes per 100kb bins versus eccDNA forming regions (averages per scaffold)
awk '{print $2}' ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage > ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage.countcolumn
paste ${SAMPLE}.genesper100kb.scaffoldaverage ${SAMPLE}.normalized.splitreadsper100kb.scaffoldaverage.countcolumn > ${SAMPLE}.SRsvsgenesperk100kbperscaffold