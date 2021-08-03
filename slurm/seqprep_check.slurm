#!/bin/bash
#SBATCH --job-name=seqprep_check
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/human/seqprep_check

mapfile="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.contignames"
genome_bwa="/global/scratch/users/pierrj/eccDNA/2018_moller/references/GRCh37.p13_bwa"
sample=SRR6315404
SAMPLE=SRR6315404

/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${SAMPLE} -O .
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e ${THREADS} -O . -t tmp ${SAMPLE}.sra


/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f ${READONE} -r ${READTWO} -1 seqprep.${sample}_R1.fastq.gz -2 seqprep.${sample}_R2.fastq.gz -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s merged.${sample}.fastq.gz
cutadapt -j ${THREADS} --nextseq-trim=20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o trimmed.seqprep.${sample}_R1.fastq -p trimmed.seqprep.${sample}_R2.fastq seqprep.${sample}_R1.fastq.gz seqprep.${sample}_R2.fastq.gz
cutadapt -j ${THREADS} --nextseq-trim=20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o trimmed.merged.${sample} merged.${sample}.fastq.gz
bwa mem -q -a -t ${THREADS} ${GENOME_DB} trimmed.seqprep.${sample}_R1.fastq trimmed.seqprep.${sample}_R2.fastq -o seqprep.trimmed.${sample}_bwamem.sam
bwa mem -q -a -t ${THREADS} ${GENOME_DB} trimmed.merged.${sample} -o merged.trimmed.${sample}_bwamem.sam
samtools view -S -b seqprep.trimmed.${sample}_bwamem.sam > seqprep.trimmed.${sample}_bwamem.bam
samtools view -S -b merged.trimmed.${sample}_bwamem.sam > merged.trimmed.${sample}_bwamem.bam
bamtools merge -in seqprep.trimmed.${sample}_bwamem.bam -in merged.trimmed.${sample}_bwamem.bam -out mergedandpe.${sample}_bwamem.bam

# sort and index for filtering
samtools sort -@ ${THREADS} mergedandpe.${sample}_bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

# filter using names of contigs in mapfile
samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam
samtools view -b -F 256 filtered.sorted.${SAMPLE}.bam > no_secondary.filtered.sorted.${SAMPLE}.bam
samtools view -b -q 1 no_secondary.filtered.sorted.${SAMPLE}.bam > uniq.filtered.sorted.${SAMPLE}.bam

# COMMENTS MISSING
samtools sort -n -@ ${THREADS} filtered.sorted.${SAMPLE}.bam > multimapped.filtered.name_sorted.${SAMPLE}.bam


# remove tmp files and unfiltered files
# tmp command should probably be more specific, or use a tmp directory
rm tmp*
rm ${SAMPLE}.mergedandpe.bwamem.bam
rm filtered.sorted.${SAMPLE}.bam
rm ${SAMPLE}.sorted.mergedandpe.bwamem.bam

/global/home/users/pierrj/git/bash/ecc_caller_mapq0.sh -m ${mapfile} \
    -s ${sample} \
    -t ${SLURM_NTASKS} \
    -b uniq.filtered.sorted.${sample}.bam \
    -q multimapped.filtered.name_sorted.${sample}.bam

/global/home/users/pierrj/git/bash/assign_confidence_nodb_nomerge.sh -m ${mapfile} \
    -s ${sample} -t ${SLURM_NTASKS} -b no_secondary.filtered.sorted.${sample}.bam -r ${sample}.confirmedsplitreads.bed