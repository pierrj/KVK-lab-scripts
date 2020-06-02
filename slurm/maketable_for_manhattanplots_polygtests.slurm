#!/bin/bash
#SBATCH --job-name=maketable_for_manhattanplots_polygtests
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/polygtrimming
table=manhattanplots
sample=IF_1A
normalize_factor=$(samtools view -c -F 4 -F 2048 sorted.MG8.mergedandpe.${sample}_bwamem.bam | awk '{print $1/1000000}' )
awk -v N=$normalize_factor -v OFS='\t' '{sum+=$3} NR%100==0 {print $1, $2}' genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam > ${table}.firstcolumn.100
awk -v N=$normalize_factor -v OFS='\t' '{sum+=$3} NR%1000==0 {print $1, $2}' genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam > ${table}.firstcolumn.1000
awk -v N=$normalize_factor '{sum+=$3} NR%100==0 {print sum/100/N; sum =0}' /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam > nopolyg.${sample}.${table}.100
awk -v N=$normalize_factor '{sum+=$3} NR%1000==0 {print sum/1000/N; sum =0}' /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam > nopolyg.${sample}.${table}.1000
awk -v N=$normalize_factor '{sum+=$3} NR%100==0 {print sum/100/N; sum =0}' genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam > yespolyg.${sample}.${table}.100
awk -v N=$normalize_factor '{sum+=$3} NR%1000==0 {print sum/1000/N; sum =0}' genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam > yespolyg.${sample}.${table}.1000

echo -e 'CHROMOSOME''\t''BASE''\t''COUNT' > ${table}.first_row
paste ${table}.firstcolumn.100 nopolyg.${sample}.${table}.100 > nopolyg.${sample}.${table}.100.table.old
paste ${table}.firstcolumn.1000 nopolyg.${sample}.${table}.1000 > nopolyg.${sample}.${table}.1000.table.old
paste ${table}.firstcolumn.100 yespolyg.${sample}.${table}.100 > yespolyg.${sample}.${table}.100.table.old
paste ${table}.firstcolumn.1000 yespolyg.${sample}.${table}.1000 > yespolyg.${sample}.${table}.1000.table.old

cat ${table}.first_row nopolyg.${sample}.${table}.100.table.old > nopolyg.${sample}.${table}.100.table
cat ${table}.first_row nopolyg.${sample}.${table}.1000.table.old > nopolyg.${sample}.${table}.1000.table
cat ${table}.first_row yespolyg.${sample}.${table}.100.table.old > yespolyg.${sample}.${table}.100.table
cat ${table}.first_row yespolyg.${sample}.${table}.1000.table.old > yespolyg.${sample}.${table}.1000.table