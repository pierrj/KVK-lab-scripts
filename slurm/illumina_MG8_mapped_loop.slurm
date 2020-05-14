#!/bin/bash
#SBATCH --job-name=doesntmatter
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
cd $SLURM_SUBMIT_DIR
cd ${sample}
/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f ${sample}_R1.fastq -r ${sample}_R2.fastq -1 tmp.seqprep.${sample}_R1.fastq.gz -2 tmp.seqprep.${sample}_R2.fastq.gz -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s tmp.merged.${sample}.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o tmp.trimmed.seqprep.${sample}_R1.fastq -p tmp.trimmed.seqprep.${sample}_R2.fastq tmp.seqprep.${sample}_R1.fastq.gz tmp.seqprep.${sample}_R2.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o tmp.trimmed.merged.${sample} tmp.merged.${sample}.fastq.gz
bwa mem -t 24 /global/scratch/users/pierrj/references/Magnaporthe_oryzae.MG8.dna.toplevel.noMT_bwa tmp.trimmed.seqprep.${sample}_R1.fastq tmp.trimmed.seqprep.${sample}_R2.fastq -o tmp.seqprep.trimmed.${sample}_bwamem.sam
bwa mem -t 24 /global/scratch/users/pierrj/references/Magnaporthe_oryzae.MG8.dna.toplevel.noMT_bwa tmp.trimmed.merged.${sample} -o tmp.merged.trimmed.${sample}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${sample}_bwamem.sam > tmp.seqprep.trimmed.${sample}_bwamem.bam
samtools view -S -b tmp.merged.trimmed.${sample}_bwamem.sam > tmp.merged.trimmed.${sample}_bwamem.bam
bamtools merge -in tmp.seqprep.trimmed.${sample}_bwamem.bam -in tmp.merged.trimmed.${sample}_bwamem.bam -out mergedandpe.${sample}_bwamem.bam
samtools sort mergedandpe.${sample}_bwamem.bam > sorted.mergedandpe.${sample}_bwamem.bam
samtools index sorted.mergedandpe.${sample}_bwamem.bam

bedtools genomecov -d -ibam sorted.mergedandpe.${sample}_bwamem.bam > genomecoverage.mergedandpe.${sample}_bwamem.bed