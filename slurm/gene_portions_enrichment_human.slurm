#!/bin/bash
#SBATCH --job-name=gene_portions_enrichment_human
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human

GENE_GFF=/global/scratch/users/pierrj/references/GRCh37.p13.gff
GENE_BEDFILE=/global/scratch/users/pierrj/references/GRCh37.p13.genes_pseudogenes
GENOME_CHROMSIZES=/global/scratch/pierrj/references/GRCh37.p13.chromsizes
SAMPLE_MAPFILE=sample_mapfile
GENOME_FILE=/global/scratch/users/pierrj/references/GRCh37.p13.fasta
copia_file=/global/scratch/pierrj/references/te_annotations/human/blank
gypsy_file=/global/scratch/pierrj/references/te_annotations/human/blank
LTR_TE_LOCATIONS=/global/scratch/pierrj/references/te_annotations/human/blank

## first make bed files

# exons

awk '$3 ~ /exon/' ${GENE_GFF} | awk -v OFS='\t' '{print $1, $4, $5, $9}' > exons

# first exon

grep "1;P" exons > first_exons

# introns

module load genometools

gt gff3 -addintrons ${GENE_GFF} | awk '$3 ~ /intron/' | awk -v OFS='\t' '{print $1, $4, $5, $9}' > introns

# 1000 bp upstream

awk '$3 ~ /gene/' ${GENE_GFF} | awk -v OFS='\t' '{if ($7 == "+")
    {print $1, $4-1000, $4}
    else if ($7 == "-")
    {print $1, $5, $5+1000}}' | awk '$2 > 0' > upstream_old
    
## to deal with cutoffs past end of genome
bedtools slop -b 0 -i upstream_old -g ${GENOME_CHROMSIZES} | awk -v OFS='\t' '{print $1, $2, $3, NR}' > upstream

# 1000 bp downstream

awk '$3 == "gene"' ${GENE_GFF} | awk -v OFS='\t' '{if ($7 == "-")
    {print $1, $4-1000, $4}
    else if ($7 == "+")
    {print $1, $5, $5+1000}}' | awk '$2 > 0' > downstream_old
    
## to deal with cutoffs past end of genome
bedtools slop -b 0 -i downstream_old -g ${GENOME_CHROMSIZES} | awk -v OFS='\t' '{print $1, $2, $3, NR}' > downstream

# genic

cp ${GENE_BEDFILE} genic

# nongenic

bedtools complement -i genic -g ${GENOME_CHROMSIZES} | awk -v OFS='\t' '{print $1, $2, $3, NR}' > intergenic

# cpg islands

module load emboss

cpgplot -sequence ${GENOME_FILE} \
        -window 100 \
        -minlen 200 \
        -minoe 0.6 \
        -minpc 50. \
        -outfeat cpg_islands.gff \
        -outfile human.plot \
        -noplot -nocg -nopc -noobsexp

awk -v OFS='\t' '{if ($1 !~ /^#/) {print $1, $4, $5, $9}}' cpg_islands.gff > cpg_islands

# first 100 bp of genes

awk -v OFS='\t' '{print $1, $2, $2+100, $4}' ${GENE_BEDFILE} > gene_first_100bp

# 100 bp starting 1000 bp upstream

awk '$3 == "gene"' ${GENE_GFF} | awk -v OFS='\t' '{if ($7 == "+")
    {print $1, $4-1000, $4-900}
    else if ($7 == "-")
    {print $1, $5+900, $5+1000}}' | awk '$2 > 0' > 100bp_starting_1000bp_upstream_old
    
## to deal with cutoffs past end of genome
bedtools slop -b 0 -i 100bp_starting_1000bp_upstream_old -g ${GENOME_CHROMSIZES} | awk -v OFS='\t' '{print $1, $2, $3, NR}' > 100bp_starting_1000bp_upstream


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
echo first_exons >> feature_mapfile
echo gene_first_100bp >> feature_mapfile
echo 100bp_starting_1000bp_upstream >> feature_mapfile

## for large eccdnas

ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human/human.large_eccdnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/

while read sample; do
cd ${sample}
    bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${copia_file} | \
        awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
        awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}' | \
        bedtools intersect -wao -a - -b ${gypsy_file} | \
        awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
        awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}'| sort -k1,1 -k2,2n | awk '$3-$2 > 400' > large_eccdnas
    sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' large_eccdnas - > large_eccdnas_splitreads
    cp large_eccdnas_splitreads ${sample}.large_eccdnas_splitreads
    realpath ${sample}.large_eccdnas_splitreads >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile_muscle

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human/

if [ -f "observed_averages_large" ]; then
    rm observed_averages_large
fi

while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn" ]; then
        rm ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
    fi
    while read ECCDNA_FILE; do
        ecc_basename=$(basename ${ECCDNA_FILE})
        #bedtools coverage -a ${FEATURE_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
        bedtools intersect -wa -c -a ${FEATURE_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
        num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
        awk -v N=$num_srs '{print $1, $2/N}' ${ecc_basename}.splitreadsperfeature > ${ecc_basename}.normalized.splitreadsperfeature
        echo ${ecc_basename}.normalized.splitreadsperfeature >> ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
    done < ${ECCDNA_MAPFILE}

    # normalize and average across technical and biological replicates as written in previous scripts
    if [ -f "${FEATURE_FILE}.normalize_table_column" ]; then
        rm ${FEATURE_FILE}.normalize_table_column
    fi
    sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
    for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${FEATURE_FILE}.normalize_table_column ; done
    paste ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn ${FEATURE_FILE}.normalize_table_column ${SAMPLE_MAPFILE} > ${FEATURE_FILE}.mapfile_for_normalize_and_average
    /global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${FEATURE_FILE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
    mv G3.normalized_binned ${FEATURE_FILE}.normalized.splitreadsperfeature

    awk -v f=${FEATURE_FILE} -v OFS='\t' '{ sum += $2} END {print f, sum/NR}' ${FEATURE_FILE}.normalized.splitreadsperfeature >> observed_averages_large
done < feature_mapfile

# shuffled

while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.permuted" ]; then
        rm ${FEATURE_FILE}.permuted
    fi
done < feature_mapfile

for i in {0..2}; do
    if [ -f "permuted_ecc_mapfile" ]; then
        rm permuted_ecc_mapfile
    fi
    while read ECCDNA_FILE; do
        bedtools shuffle -i ${ECCDNA_FILE} -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > shuffled.${ecc_basename}
        echo shuffled.${ecc_basename} >> permuted_ecc_mapfile
    done < ${ECCDNA_MAPFILE}
    while read FEATURE_FILE; do
        if [ -f "${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn" ]; then
            rm ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
        fi
        while read ECCDNA_FILE; do
            ecc_basename=$(basename ${ECCDNA_FILE})
            #bedtools coverage -a ${FEATURE_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
            bedtools intersect -wa -c -a ${FEATURE_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
            num_srs=$(wc -l ${ecc_basename} | awk '{print $1/100000}')
            awk -v N=$num_srs '{print $1, $2/N}' ${ecc_basename}.splitreadsperfeature > ${ecc_basename}.normalized.splitreadsperfeature
            echo ${ecc_basename}.normalized.splitreadsperfeature >> ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
        done < permuted_ecc_mapfile

        # normalize and average across technical and biological replicates as written in previous scripts
        if [ -f "${FEATURE_FILE}.normalize_table_column" ]; then
            rm ${FEATURE_FILE}.normalize_table_column
        fi
        sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
        for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${FEATURE_FILE}.normalize_table_column ; done
        paste ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn ${FEATURE_FILE}.normalize_table_column ${SAMPLE_MAPFILE} > ${FEATURE_FILE}.mapfile_for_normalize_and_average
        /global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${FEATURE_FILE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
        mv G3.normalized_binned ${FEATURE_FILE}.normalized.splitreadsperfeature

        awk '{ sum += $2} END {print sum/NR}' ${FEATURE_FILE}.normalized.splitreadsperfeature >> ${FEATURE_FILE}.permuted
    done < feature_mapfile
done

if [ -f "expected_averages_large" ]; then
    rm expected_averages_large
fi

while read FEATURE_FILE; do
    awk -v f=${FEATURE_FILE} -v OFS='\t' '{ sum += $1} END {print f, sum/NR}' ${FEATURE_FILE}.permuted >> expected_averages_large
done < feature_mapfile

## for microdnas

ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human/human.micro_dnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/2018_moller/full_run/

while read sample; do
cd ${sample}
    bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b ${copia_file} | \
        awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
        awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}' | \
        bedtools intersect -wao -a - -b ${gypsy_file} | \
        awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
        awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}'| sort -k1,1 -k2,2n | awk '$3-$2 <= 400' > ${sample}.micro_dnas
    sed 's/[[:space:]]*$//' ${sample}.ecc_caller_out.splitreads.bed | awk 'NR==FNR{a[$0]++;next}a[$0]' micro_dnas - > micro_dnas_splitreads
    cp micro_dnas_splitreads ${sample}.micro_dnas_splitreads
    realpath ${sample}.micro_dnas_splitreads >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile_muscle

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/human/

if [ -f "observed_averages_micro" ]; then
    rm observed_averages_micro
fi

while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn" ]; then
        rm ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
    fi
    while read ECCDNA_FILE; do
        ecc_basename=$(basename ${ECCDNA_FILE})
        #bedtools coverage -a ${FEATURE_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
        bedtools intersect -wa -c -a ${FEATURE_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
        num_srs=$(wc -l ${ECCDNA_FILE} | awk '{print $1/100000}')
        awk -v N=$num_srs '{print $1, $2/N}' ${ecc_basename}.splitreadsperfeature > ${ecc_basename}.normalized.splitreadsperfeature
        echo ${ecc_basename}.normalized.splitreadsperfeature >> ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
    done < ${ECCDNA_MAPFILE}

    # normalize and average across technical and biological replicates as written in previous scripts
    if [ -f "${FEATURE_FILE}.normalize_table_column" ]; then
        rm ${FEATURE_FILE}.normalize_table_column
    fi
    sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
    for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${FEATURE_FILE}.normalize_table_column ; done
    paste ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn ${FEATURE_FILE}.normalize_table_column ${SAMPLE_MAPFILE} > ${FEATURE_FILE}.mapfile_for_normalize_and_average
    /global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${FEATURE_FILE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
    mv G3.normalized_binned ${FEATURE_FILE}.normalized.splitreadsperfeature

    awk -v f=${FEATURE_FILE} -v OFS='\t' '{ sum += $2} END {print f, sum/NR}' ${FEATURE_FILE}.normalized.splitreadsperfeature >> observed_averages_micro
done < feature_mapfile

# shuffled

while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.permuted" ]; then
        rm ${FEATURE_FILE}.permuted
    fi
done < feature_mapfile

for i in {0..2}; do
    if [ -f "permuted_ecc_mapfile" ]; then
        rm permuted_ecc_mapfile
    fi
    while read ECCDNA_FILE; do
        bedtools shuffle -i ${ECCDNA_FILE} -g ${GENOME_CHROMSIZES} -excl ${LTR_TE_LOCATIONS} > shuffled.${ecc_basename}
        echo shuffled.${ecc_basename} >> permuted_ecc_mapfile
    done < ${ECCDNA_MAPFILE}
    while read FEATURE_FILE; do
        if [ -f "${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn" ]; then
            rm ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
        fi
        while read ECCDNA_FILE; do
            ecc_basename=$(basename ${ECCDNA_FILE})
            #bedtools coverage -a ${FEATURE_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
            bedtools intersect -wa -c -a ${FEATURE_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
            num_srs=$(wc -l ${ecc_basename} | awk '{print $1/100000}')
            awk -v N=$num_srs '{print $1, $2/N}' ${ecc_basename}.splitreadsperfeature > ${ecc_basename}.normalized.splitreadsperfeature
            echo ${ecc_basename}.normalized.splitreadsperfeature >> ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
        done < permuted_ecc_mapfile

        # normalize and average across technical and biological replicates as written in previous scripts
        if [ -f "${FEATURE_FILE}.normalize_table_column" ]; then
            rm ${FEATURE_FILE}.normalize_table_column
        fi
        sample_count=$(wc -l ${SAMPLE_MAPFILE} | awk '{print $1+1}')
        for (( i = 1 ; i < ${sample_count}; i++)); do echo 1 >> ${FEATURE_FILE}.normalize_table_column ; done
        paste ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn ${FEATURE_FILE}.normalize_table_column ${SAMPLE_MAPFILE} > ${FEATURE_FILE}.mapfile_for_normalize_and_average
        /global/home/users/pierrj/git/bash/normalize_and_average.sh -m ${FEATURE_FILE}.mapfile_for_normalize_and_average -f 1 -b 1 -c 2 -n n
        mv G3.normalized_binned ${FEATURE_FILE}.normalized.splitreadsperfeature

        awk '{ sum += $2} END {print sum/NR}' ${FEATURE_FILE}.normalized.splitreadsperfeature >> ${FEATURE_FILE}.permuted
    done < feature_mapfile
done

if [ -f "expected_averages_micro" ]; then
    rm expected_averages_micro
fi

while read FEATURE_FILE; do
    awk -v f=${FEATURE_FILE} -v OFS='\t' '{ sum += $1} END {print f, sum/NR}' ${FEATURE_FILE}.permuted >> expected_averages_micro
done < feature_mapfile