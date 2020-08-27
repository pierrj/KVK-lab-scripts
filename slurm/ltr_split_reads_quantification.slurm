#!/bin/bash
#SBATCH --job-name=ltr_split_reads_quantification
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ltr_splitreads_tests/


sample=CT_1A_pe_noseqprep
bamfile=/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/filtered.sorted.CT_1A.bam

/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -s ${sample} -b ${bamfile}

sample=CT_1A_se_noseqprep
bamfile=/global/scratch/users/pierrj/eccDNA/pipeline_tests/single_end_test/filtered.sorted.CT_1A.bam

/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -s ${sample} -b ${bamfile}

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ltr_splitreads_tests/seqprep_singleend
SAMPLE="CT_1A"
READONE="/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/CT_1A_R1.fastq"
READTWO="/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/CT_1A_R2.fastq"
GENOME_DB="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_with_spikeins_tweaked_bwa"
THREADS="24"
MAPFILE="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"

/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f ${READONE} -r ${READTWO} -1 tmp.seqprep.${SAMPLE}_R1.fastq.gz -2 tmp.seqprep.${SAMPLE}_R2.fastq.gz -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s tmp.merged.${SAMPLE}.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq tmp.seqprep.${SAMPLE}_R1.fastq.gz tmp.seqprep.${SAMPLE}_R2.fastq.gz
cat tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq > tmp.trimmed.seqprep.${SAMPLE}.fastq
bwa mem -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
bwa mem -t ${THREADS} ${GENOME_DB} tmp.trimmed.merged.${SAMPLE}.fastq -o tmp.merged.trimmed.${SAMPLE}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > tmp.seqprep.trimmed.${SAMPLE}_bwamem.bam
samtools view -S -b tmp.merged.trimmed.${SAMPLE}_bwamem.sam > tmp.merged.trimmed.${SAMPLE}_bwamem.bam
bamtools merge -in tmp.seqprep.trimmed.${SAMPLE}_bwamem.bam -in tmp.merged.trimmed.${SAMPLE}_bwamem.bam -out ${SAMPLE}.mergedandpe.bwamem.bam
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam

cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/ltr_splitreads_tests/

sample=CT_1A_se_yesseqprep
bamfile=/global/scratch/users/pierrj/eccDNA/pipeline_tests/ltr_splitreads_tests/seqprep_singleend/filtered.sorted.${SAMPLE}.bam

/global/home/users/pierrj/git/bash/get_splitreads_for_LTR.sh -s ${sample} -b ${bamfile}