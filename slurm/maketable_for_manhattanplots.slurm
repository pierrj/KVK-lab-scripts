#!/bin/bash
#SBATCH --job-name=maketable_for_manhattanplots
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina
table=manhattanplots
while read sample; 
do
    cd ${sample}
    normalize_factor=$(samtools view -c -F 4 -F 2048 sorted.MG8.mergedandpe.${sample}_bwamem.bam | awk '{print $1/1000000}' )
    awk -v N=$normalize_factor '{sum+=$3} NR%100==0 {print sum/100/N; sum =0}' genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam > /global/scratch/users/pierrj/eccDNA/pipeline_tests/manhattanplots/${sample}.${table}
    cd ..
done < mapfile
sample=CT_1A
cd CT_1A
normalize_factor=$(samtools view -c -F 4 -F 2048 sorted.MG8.mergedandpe.${sample}_bwamem.bam | awk '{print $1/1000000}' )
awk -v N=$normalize_factor -v OFS='\t' '{sum+=$3} NR%100==0 {print $1, $2}' /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam > /global/scratch/users/pierrj/eccDNA/pipeline_tests/manhattanplots/${table}.firstcolumn
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/manhattanplots/
if [ ! -d "techrepsmerged_${table}" ]; then
    mkdir techrepsmerged_${table}
fi
ls *.${table} | awk '{print substr($0,0,4)}' | sort | uniq > techrepsmerged_${table}/mapfile
while read sample; 
do
    paste $(find . -maxdepth 1 -name "${sample}*.${table}" | xargs -r ls -1 | cut -c 3- | tr "\n" " ") | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' > techrepsmerged_${table}/${sample}.${table}
done < techrepsmerged_${table}/mapfile
cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/manhattanplots/${table}.firstcolumn /global/scratch/users/pierrj/eccDNA/pipeline_tests/manhattanplots/techrepsmerged_${table}/treatmentsmerged_${table}/${table}.firstcolumn
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
cp /global/scratch/users/pierrj/eccDNA/pipeline_tests/manhattanplots/${table}.firstcolumn /global/scratch/users/pierrj/eccDNA/pipeline_tests/manhattanplots/techrepsmerged_${table}/treatmentsmerged_${table}/${table}.firstcolumn
while read sample; 
do
    paste ${table}.firstcolumn ${sample}.${table} > ${sample}.${table}.table.old
    echo -e 'CHROMOSOME''\t''BASE''\t''COUNT' > ${table}.first_row
    cat ${table}.first_row ${sample}.${table}.table.old > ${sample}.${table}.table
done < mapfile