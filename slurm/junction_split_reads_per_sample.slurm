#!/bin/bash
#SBATCH --job-name=junction_split_reads_per_sample
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina

output_table=/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/srs_per_reads

if [ -f "${output_table}" ]; then
    rm ${output_table}
fi

while read sample; do
    cd ${sample}
        read_count=$(samtools view -c no_secondary.filtered.sorted.${sample}.bam)
        sr_count=$(wc -l ${sample}.ecc_caller_out.splitreads.bed | awk '{print $1}')
        echo -e ${sample}'\t'${read_count}'\t'${sr_count} >> /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/srs_per_reads
    cd ..
done < mapfile

awk '{print $3/$2}' ${output_table} | awk '{SUM+=$1} END {print SUM/NR}'