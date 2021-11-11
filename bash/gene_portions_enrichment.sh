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
#SOFTWARE.
#!/bin/bash
while getopts o:g:s:d:f:c:y:n: option
do
case "${option}"
in
o) ORGANISM_NAME=${OPTARG};;
g) GENE_GFF=${OPTARG};;
s) SAMPLE_MAPFILE=${OPTARG};;
d) ECC_DIR=${OPTARG};;
f) GENOME_FILE=${OPTARG};;
c) COPIA_FILE=${OPTARG};;
y) GYPSY_FILE=${OPTARG};;
n) CONTIGNAMES_FILE=${OPTARG};;
esac
done

WORKING_DIR=$(pwd)

## USAGE ##
# this script calculates which regions of the genome are enriched for micro DNAs and large eccDNAs
# -o name of the organism
# -g gff file of gene annotations
# -s mapfile containing list of eccDNA sample names
# -d directory where eccDNA data is locatied
# -f genome fasta file
# -c bed file of copia element locations
# -y bed file of gypsy element locations
# -n list of scaffolds of interest within genome file


# first off filter all files to scaffolds of interest, the ones use for eccdna calling
samtools faidx ${GENOME_FILE} $(cat ${CONTIGNAMES_FILE} | tr "\n" " ") > ${GENOME_FILE}.filtered
GENOME_FILE=$(realpath ${GENOME_FILE}.filtered)

samtools faidx ${GENOME_FILE}
cut -f1,2 ${GENOME_FILE}.fai > ${GENOME_FILE}.chromsizes
GENOME_CHROMSIZES=$(realpath ${GENOME_FILE}.chromsizes)

awk -v OFS='\t' '{print $1, 1, $2}' ${GENOME_CHROMSIZES} > ${GENOME_CHROMSIZES}.forintersect
GENOME_CHROMSIZES_FOR_INTERSECT=$(realpath ${GENOME_CHROMSIZES}.forintersect)

bedtools intersect -header -u -a ${GENE_GFF} -b ${GENOME_CHROMSIZES_FOR_INTERSECT} > ${GENE_GFF}.filtered
GENE_GFF=$(realpath ${GENE_GFF}.filtered)

bedtools intersect -header -u -a ${COPIA_FILE} -b ${GENOME_CHROMSIZES_FOR_INTERSECT} > ${COPIA_FILE}.filtered
COPIA_FILE=$(realpath ${COPIA_FILE}.filtered)

bedtools intersect -header -u -a ${GYPSY_FILE} -b ${GENOME_CHROMSIZES_FOR_INTERSECT} > ${GYPSY_FILE}.filtered
GYPSY_FILE=$(realpath ${GYPSY_FILE}.filtered)


# make a gene bedfile
basename_gff=$(basename ${GENE_GFF})
awk -v OFS='\t' '{if ($3 ~ /gene/) {print $1, $4, $5, $9}}' ${GENE_GFF} | sort -k1,1 -k2,2n > ${basename_gff}.bed
GENE_BEDFILE=$(realpath ${basename_gff}.bed)

# all te locations
cat ${COPIA_FILE} ${GYPSY_FILE} > ${ORGANISM}.ltr_te_locs
LTR_TE_LOCS=$(realpath ${ORGANISM}.ltr_te_locs)

# get lengths of scaffolds of interest
samtools faidx ${GENOME_FILE}
cut -f1,2 ${GENOME_FILE}.fai > ${GENOME_FILE}.chromsizes
GENOME_CHROMSIZES=$(realpath ${GENOME_FILE}.chromsizes)

# fix and standardize gffs
agat_convert_sp_gff2gtf.pl --gff ${GENE_GFF} -o ${ORGANISM}.gtf > agat_output 2>&1
GENE_GTF=$(realpath ${ORGANISM}.gtf)



## make regions of interest files

# exons

awk '$3 ~ /exon/' ${GENE_GFF} | awk -v OFS='\t' '{print $1, $4, $5, $9}' > exons

# introns

agat_sp_add_introns.pl --gff ${GENE_GFF} -o ${ORGANISM}.tmpintrons > agat_output_intron 2>&1

## for some reason human gff file has a weird intron so remove it here
awk '$3 ~ /intron/' ${ORGANISM}.tmpintrons | awk -v OFS='\t' '{ if ($5 > $4) {print $1, $4, $5, $9}}' > introns

# 2000 bp upstream

awk '$3 ~ /gene/' ${GENE_GFF} | awk -v OFS='\t' '{if ($7 == "+")
    {print $1, $4-2000, $4}
    else if ($7 == "-")
    {print $1, $5, $5+2000}}' | awk '$2 > 0' > upstream_old
    
## to deal with cutoffs past end of genome
bedtools slop -b 0 -i upstream_old -g ${GENOME_CHROMSIZES} | awk -v OFS='\t' '{print $1, $2, $3, NR}' > upstream

# 1000 bp downstream

awk '$3 == "gene"' ${GENE_GFF} | awk -v OFS='\t' '{if ($7 == "-")
    {print $1, $4-2000, $4}
    else if ($7 == "+")
    {print $1, $5, $5+2000}}' | awk '$2 > 0' > downstream_old
    
## to deal with cutoffs past end of genome
bedtools slop -b 0 -i downstream_old -g ${GENOME_CHROMSIZES} | awk -v OFS='\t' '{print $1, $2, $3, NR}' > downstream

# genic

cp ${GENE_BEDFILE} genic

# nongenic

bedtools complement -i genic -g ${GENOME_CHROMSIZES} | awk -v OFS='\t' '{print $1, $2, $3, NR}' > intergenic

# cpg islands

cpgplot -sequence ${GENOME_FILE} \
        -window 100 \
        -minlen 200 \
        -minoe 0.6 \
        -minpc 50. \
        -outfeat cpg_islands.gff \
        -outfile human.plot \
        -noplot -nocg -nopc -noobsexp

awk -v OFS='\t' '{if ($1 !~ /^#/) {print $1, $4, $5, $9}}' cpg_islands.gff > cpg_islands

# five_prime_utr

awk -v OFS='\t' '{ if ($3 ~ /five_prime_utr/) {print $1, $4, $5, NR}}' ${GENE_GTF} > five_prime_utr

# three_prime_utr 
awk -v OFS='\t' '{ if ($3 ~ /three_prime_utr/) {print $1, $4, $5, NR}}' ${GENE_GTF} > three_prime_utr

# make mapfile to loop through
if [ -f "feature_mapfile" ]; then
    rm feature_mapfile
fi
echo exons >> feature_mapfile
echo introns >> feature_mapfile
echo upstream >> feature_mapfile
echo downstream >> feature_mapfile
echo genic >> feature_mapfile
echo intergenic >> feature_mapfile
echo cpg_islands >> feature_mapfile
echo five_prime_utr >> feature_mapfile
echo three_prime_utr >> feature_mapfile


if [ -f "large_eccdna_mapfile" ]; then
    rm large_eccdna_mapfile
fi

touch large_eccdna_mapfile
ECCDNA_MAPFILE=$(realpath large_eccdna_mapfile)

cd ${ECC_DIR}

# get large eccdna files
while read sample; do
cd ${sample}
    bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${COPIA_FILE} | \
        awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
        awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $1, $2, $3}}' | \
        bedtools intersect -wao -a - -b ${GYPSY_FILE} | \
        awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
        awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $1, $2, $3}}'| sort -k1,1 -k2,2n | awk '$3-$2 > 400' > large_eccdnas
    sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' large_eccdnas - > large_eccdnas_splitreads
    cp large_eccdnas_splitreads ${sample}.large_eccdnas_splitreads
    realpath ${sample}.large_eccdnas_splitreads >> ${ECCDNA_MAPFILE}
cd ..
done < ${SAMPLE_MAPFILE}

cd ${WORKING_DIR}

if [ -f "observed_averages_large" ]; then
    rm observed_averages_large
fi
# for each find the percentage of large eccdnas with any overlap
while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.percent_per_sample" ]; then
        rm ${FEATURE_FILE}.percent_per_sample
    fi
    while read ECCDNA_FILE; do
        num_intersect=$(bedtools intersect -u -a ${ECCDNA_FILE} -b ${FEATURE_FILE} | wc -l )
        num_total=$(wc -l ${ECCDNA_FILE} | awk '{print $1}')
        percentage=$(awk -v var1=$num_intersect -v var2=$num_total 'BEGIN { print  ( var1 / var2 ) }')
        echo -e $(basename ${ECCDNA_FILE})'\t'${percentage} >> ${FEATURE_FILE}.percent_per_sample
    done < ${ECCDNA_MAPFILE}
    awk -v f=${FEATURE_FILE} -v OFS='\t' '{sum+=$2} END {print f, sum/NR}' ${FEATURE_FILE}.percent_per_sample >> observed_averages_large
done < feature_mapfile

while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.permuted" ]; then
        rm ${FEATURE_FILE}.permuted
    fi
done < feature_mapfile

# generate simulated eccdnas
for i in {0..9}; do
    if [ -f "permuted_ecc_mapfile" ]; then
        rm permuted_ecc_mapfile
    fi
    while read ECCDNA_FILE; do
        ecc_basename=$(basename ${ECCDNA_FILE})
        # shuffle eccdnas, excluding ltr retrotranposon locations
        bedtools shuffle -i ${ECCDNA_FILE} -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > shuffled.${ecc_basename}
        echo shuffled.${ecc_basename} >> permuted_ecc_mapfile
    done < ${ECCDNA_MAPFILE}
    # same as before but with shuffled eccdnas
    while read FEATURE_FILE; do
        if [ -f "${FEATURE_FILE}.percent_per_sample" ]; then
            rm ${FEATURE_FILE}.percent_per_sample
        fi
        while read ECCDNA_FILE; do
            num_intersect=$(bedtools intersect -u -a ${ECCDNA_FILE} -b ${FEATURE_FILE} | wc -l )
            num_total=$(wc -l ${ECCDNA_FILE} | awk '{print $1}')
            percentage=$(awk -v var1=$num_intersect -v var2=$num_total 'BEGIN { print  ( var1 / var2 ) }')
            echo -e $(basename ${ECCDNA_FILE})'\t'${percentage} >> ${FEATURE_FILE}.percent_per_sample
        done < permuted_ecc_mapfile
        awk -v f=${FEATURE_FILE} -v OFS='\t' '{sum+=$2} END {print sum/NR}' ${FEATURE_FILE}.percent_per_sample >> ${FEATURE_FILE}.permuted
    done < feature_mapfile
done

if [ -f "expected_averages_large" ]; then
    rm expected_averages_large
fi

# final expected output
while read FEATURE_FILE; do
    awk -v f=${FEATURE_FILE} -v OFS='\t' '{ sum += $1} END {print f, sum/NR}' ${FEATURE_FILE}.permuted >> expected_averages_large
done < feature_mapfile

if [ -f "micro_dna_mapfile" ]; then
    rm micro_dna_mapfile
fi

# same as above but for microdnas
touch micro_dna_mapfile
ECCDNA_MAPFILE=$(realpath micro_dna_mapfile)

cd ${ECC_DIR}

# get microdnas files
while read sample; do
cd ${sample}
    bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${COPIA_FILE} | \
        awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
        awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $1, $2, $3}}' | \
        bedtools intersect -wao -a - -b ${GYPSY_FILE} | \
        awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
        awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $1, $2, $3}}'| sort -k1,1 -k2,2n | awk '$3-$2 <= 400' > micro_dnas
    sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' micro_dnas - > micro_dnas_splitreads
    cp micro_dnas_splitreads ${sample}.micro_dnas_splitreads
    realpath ${sample}.micro_dnas_splitreads >> ${ECCDNA_MAPFILE}
cd ..
done < ${SAMPLE_MAPFILE}

cd ${WORKING_DIR}

if [ -f "observed_averages_micro" ]; then
    rm observed_averages_micro
fi

# for each find the percentage of micro dnas with any overlap
while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.percent_per_sample" ]; then
        rm ${FEATURE_FILE}.percent_per_sample
    fi
    while read ECCDNA_FILE; do
        num_intersect=$(bedtools intersect -u -a ${ECCDNA_FILE} -b ${FEATURE_FILE} | wc -l )
        num_total=$(wc -l ${ECCDNA_FILE} | awk '{print $1}')
        percentage=$(awk -v var1=$num_intersect -v var2=$num_total 'BEGIN { print  ( var1 / var2 ) }')
        echo -e $(basename ${ECCDNA_FILE})'\t'${percentage} >> ${FEATURE_FILE}.percent_per_sample
    done < ${ECCDNA_MAPFILE}
    awk -v f=${FEATURE_FILE} -v OFS='\t' '{sum+=$2} END {print f, sum/NR}' ${FEATURE_FILE}.percent_per_sample >> observed_averages_micro
done < feature_mapfile

while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.permuted" ]; then
        rm ${FEATURE_FILE}.permuted
    fi
done < feature_mapfile

# generate simulated microdnas
for i in {0..9}; do
    if [ -f "permuted_ecc_mapfile" ]; then
        rm permuted_ecc_mapfile
    fi
    while read ECCDNA_FILE; do
        ecc_basename=$(basename ${ECCDNA_FILE})
        # shuffle eccdnas, excluding ltr retrotranposon locations
        bedtools shuffle -i ${ECCDNA_FILE} -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > shuffled.${ecc_basename}
        echo shuffled.${ecc_basename} >> permuted_ecc_mapfile
    done < ${ECCDNA_MAPFILE}
    # same as before but with shuffled eccdnas
    while read FEATURE_FILE; do
        if [ -f "${FEATURE_FILE}.percent_per_sample" ]; then
            rm ${FEATURE_FILE}.percent_per_sample
        fi
        while read ECCDNA_FILE; do
            num_intersect=$(bedtools intersect -u -a ${ECCDNA_FILE} -b ${FEATURE_FILE} | wc -l )
            num_total=$(wc -l ${ECCDNA_FILE} | awk '{print $1}')
            percentage=$(awk -v var1=$num_intersect -v var2=$num_total 'BEGIN { print  ( var1 / var2 ) }')
            echo -e $(basename ${ECCDNA_FILE})'\t'${percentage} >> ${FEATURE_FILE}.percent_per_sample
        done < permuted_ecc_mapfile
        awk -v f=${FEATURE_FILE} -v OFS='\t' '{sum+=$2} END {print sum/NR}' ${FEATURE_FILE}.percent_per_sample >> ${FEATURE_FILE}.permuted
    done < feature_mapfile
done

#final expected count
if [ -f "expected_averages_micro" ]; then
    rm expected_averages_micro
fi

while read FEATURE_FILE; do
    awk -v f=${FEATURE_FILE} -v OFS='\t' '{ sum += $1} END {print f, sum/NR}' ${FEATURE_FILE}.permuted >> expected_averages_micro
done < feature_mapfile