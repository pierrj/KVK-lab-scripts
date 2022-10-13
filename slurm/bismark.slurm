#!/bin/bash
#SBATCH --job-name=bismark
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/bismark/moryzae_test

SAMPLE=SRR653493

module load bowtie2
module load samtools

# /global/scratch/users/pierrj/bismark/Bismark-0.24.0/bismark_genome_preparation \
#     --verbose --bowtie2 --path_to_aligner /global/home/groups/consultsw/sl-7.x86_64/modules/bowtie2/2.3.4.1 guy11_genome_baoetal2017/


/global/scratch/users/pierrj/bismark/Bismark-0.24.0/bismark \
        --parallel 6 \
        --temp_dir tmp \
        --output_dir ${SAMPLE}_bismark_out \
        -p 4 \
        -1 ${SAMPLE}.sra_1.fastq -2 ${SAMPLE}.sra_2.fastq \
        guy11_genome_baoetal2017/

# /global/scratch/users/pierrj/bismark/Bismark-0.24.0/deduplicate_bismark -p \
# 		--output_dir  ${SAMPLE}_bismark_out \
# 		/global/scratch/users/chandlersutherland/e14/bismark/bam_files/${f}_1_val_1_bismark_bt2_pe.bam

# /global/scratch/users/pierrj/bismark/Bismark-0.24.0//bismark_methylation_extractor -p \
#     --output ${SAMPLE}_bismark_out \
#     --ignore_r2 2 \
#     --comprehensive \
#     --parallel $SLURM_NTASKS \
#     /global/scratch/users/chandlersutherland/e14/bismark/deduplicate_bismark/${f}_1_val_1_bismark_bt2_pe.deduplicated.bam