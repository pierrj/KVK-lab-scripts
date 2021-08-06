#!/bin/bash
#SBATCH --job-name=gene_portions_enrichment
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions


GENE_GFF=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3
GENE_BEDFILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.justgenes.bed
GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes
SAMPLE_MAPFILE=/global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison_wrapper/anyoverlap_normalization/G3.sample_mapfile

## first make bed files

# exons

grep exon ${GENE_GFF} | awk -v OFS='\t' '{print $1, $4, $5, $9}' > exons

# introns

module load genometools

gt gff3 -addintrons guy11_fungap_out_12_28_20.gff3 | grep intron | awk -v OFS='\t' '{print $1, $4, $5, $9}' > introns

# 1000 bp upstream

awk '$3 == "gene"' ${GENE_GFF} | awk -v OFS='\t' '{if ($7 == "+")
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

cpgplot -sequence guy11_genome_baoetal2017.fasta \
        -window 100 \
        -minlen 200 \
        -minoe 0.6 \
        -minpc 50. \
        -outfeat cpg_islands.gff \
        -outfile guy11.plot \
        -noplot -nocg -nopc -noobsexp

awk -v OFS='\t' '{if ($1 !~ /^#/) {print $1, $4, $5, $9}}' cpg_islands.gff > cpg_islands


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


## for large eccdnas

ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/G3.large_eccdnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do
cd ${sample}
    cp large_eccdnas_splitreads ${sample}.large_eccdnas_splitreads
    realpath ${sample}.large_eccdnas_splitreads >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions

if [ -f "observed_averages_large" ]; then
    rm observed_averages_large
fi

while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn" ]; then
        rm ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
    fi
    while read ECCDNA_FILE; do
        ecc_basename=$(basename ${ECCDNA_FILE})
        bedtools intersect -F 0.5 -wa -c -a ${FEATURE_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
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

if [ -f "expected_averages_large" ]; then
    rm expected_averages_large
fi

while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.permuted" ]; then
        rm ${FEATURE_FILE}.permuted
    fi
    for i in {0..99}; do
        if [ -f "${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn" ]; then
            rm ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
        fi
        while read ECCDNA_FILE; do
            ecc_basename=$(basename ${ECCDNA_FILE})
            bedtools shuffle -i ${ECCDNA_FILE} -g ${GENOME_CHROMSIZES}> shuffled.${ecc_basename}
            bedtools intersect -F 0.5 -wa -c -a ${FEATURE_FILE} -b shuffled.${ecc_basename} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
            num_srs=$(wc -l shuffled.${ecc_basename} | awk '{print $1/100000}')
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
        awk '{ sum += $2} END {print sum/NR}' ${FEATURE_FILE}.normalized.splitreadsperfeature >> ${FEATURE_FILE}.permuted
    done

    awk -v f=${FEATURE_FILE} -v OFS='\t' '{ sum += $1} END {print sum/NR}' ${FEATURE_FILE}.permuted >> expected_averages_large
done < feature_mapfile



## for microdnas


ECCDNA_MAPFILE=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions/G3.micro_dnas_mapfile

if [ -f "${ECCDNA_MAPFILE}" ]; then
    rm ${ECCDNA_MAPFILE}
fi

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

while read sample; do
cd ${sample}
    cp micro_dnas_splitreads ${sample}.micro_dnas_splitreads
    realpath ${sample}.micro_dnas_splitreads >> ${ECCDNA_MAPFILE}
cd ..
done < mapfile

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/ecc_characteristics/gene_portions

if [ -f "observed_averages_micro" ]; then
    rm observed_averages_micro
fi

while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn" ]; then
        rm ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
    fi
    while read ECCDNA_FILE; do
        ecc_basename=$(basename ${ECCDNA_FILE})
        bedtools intersect -F 0.5 -wa -c -a ${FEATURE_FILE} -b ${ECCDNA_FILE} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
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

if [ -f "expected_averages_micro" ]; then
    rm expected_averages_micro
fi

while read FEATURE_FILE; do
    if [ -f "${FEATURE_FILE}.permuted" ]; then
        rm ${FEATURE_FILE}.permuted
    fi
    for i in {0..99}; do
        if [ -f "${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn" ]; then
            rm ${FEATURE_FILE}.mapfile_for_normalize_and_average_filecolumn
        fi
        while read ECCDNA_FILE; do
            ecc_basename=$(basename ${ECCDNA_FILE})
            bedtools shuffle -i ${ECCDNA_FILE} -g ${GENOME_CHROMSIZES}> shuffled.${ecc_basename}
            bedtools intersect -F 0.5 -wa -c -a ${FEATURE_FILE} -b shuffled.${ecc_basename} | awk -v OFS='\t' '{print $4, $5}' > ${ecc_basename}.splitreadsperfeature
            num_srs=$(wc -l shuffled.${ecc_basename} | awk '{print $1/100000}')
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
        awk '{ sum += $2} END {print sum/NR}' ${FEATURE_FILE}.normalized.splitreadsperfeature >> ${FEATURE_FILE}.permuted
    done

    awk -v f=${FEATURE_FILE} -v OFS='\t' '{ sum += $1} END {print f, sum/NR}' ${FEATURE_FILE}.permuted >> expected_averages_micro
done < feature_mapfile