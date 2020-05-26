#!/bin/bash
#SBATCH --job-name=SRcountspergene_KROJ
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina
table=srcounts
while read sample;
do
    cd ${sample}
    bedtools intersect -f 1 -wa -c -a /global/scratch/users/pierrj/references/GUY11_PacBio_merge_KROJ_renamed.genes.gff3 -b /global/scratch/users/pierrj/eccDNA/pipeline_tests/merging/parallel.confirmed.${sample}.bed | awk -v OFS='\t' '{print $9, $10}' > /global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison/genes.parallel.confirmed.${sample}.bed
    normalize_factor=$(wc -l /global/scratch/users/pierrj/eccDNA/pipeline_tests/merging/parallel.confirmed.${sample}.bed | awk '{print $1/100000}')
    awk -v N=$normalize_factor '{print $2/N}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison/genes.parallel.confirmed.${sample}.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison/${sample}.${table}
    cd ..
done < mapfile
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison
if [ ! -d "techrepsmerged_${table}" ]; then
    mkdir techrepsmerged_${table}
fi
ls *.${table} | awk '{print substr($0,0,4)}' | sort | uniq > techrepsmerged_${table}/mapfile
while read sample; 
do
    paste $(find . -maxdepth 1 -name "${sample}*.${table}" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > techrepsmerged_${table}/${sample}.${table}
done < techrepsmerged_${table}/mapfile
cd techrepsmerged_${table}
if [ ! -d "treatmentsmerged_${table}" ]; then
    mkdir treatmentsmerged_${table}
fi
ls *.${table} | awk '{print substr($0,0,2)}' | sort | uniq > treatmentsmerged_${table}/mapfile
while read sample; 
do
    paste $(find . -maxdepth 1 -name "${sample}*.${table}" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > treatmentsmerged_${table}/${sample}.${table}
done < treatmentsmerged_${table}/mapfile
cd treatmentsmerged_${table}
sample=G3_1A
normalize_factor=$(wc -l /global/scratch/users/pierrj/eccDNA/pipeline_tests/merging/parallel.confirmed.${sample}.bed | awk '{print $1/100000}')
awk -v N=$normalize_factor '{print $1}' /global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison/genes.parallel.confirmed.${sample}.bed > /global/scratch/users/pierrj/eccDNA/pipeline_tests/rnaseq_comparison/techrepsmerged_${table}/treatmentsmerged_${table}/${table}.firstcolumn
while read sample; 
do
    paste ${table}.firstcolumn ${sample}.${table} > ${sample}.${table}.table.old
    echo -e 'CHROMOSOME''\t''BASE''\t''COUNT' > ${table}.first_row
    cat ${table}.first_row ${sample}.${table}.table.old > ${sample}.${table}.table
done < mapfile