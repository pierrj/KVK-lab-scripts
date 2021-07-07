#!/bin/bash
#SBATCH --job-name=wgs_false_positives
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/wgs_data\

if [ -f "wgs_false_positives_count" ]; then
    rm wgs_false_positives_count
fi

while read sample; do 
    echo ${sample} >> wgs_false_positives_count
    wc -l ${sample}/${sample}.sra_1.fastq | awk '{print $1*2}' >> wgs_false_positives_count
    awk '$5 != "lowq"' ${sample}/${sample}.ecc_caller_out_details.txt | wc -l >> wgs_false_positives_count
    wc -l ${sample}/${sample}.ecc_caller_out.splitreads.bed >> wgs_false_positives_count
done < mapfile