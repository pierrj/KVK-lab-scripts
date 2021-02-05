#!/bin/bash
while getopts g:1:2:s:t:m:q:w: option
do
case "${option}"
in
g) GENOME_DB=${OPTARG};;
1) READONE=${OPTARG};;
2) READTWO=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
m) MAPFILE=${OPTARG};;
q) READONE_LINK=${OPTARG};;
w) READTWO_LINK=${OPTARG};;
esac
done

curl -O ${READONE_LINK}
curl -O ${READTWO_LINK}
cutadapt -j ${THREADS} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq ${READONE} ${READTWO}

# map reads
bwa mem -q -a -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam

# generate bam file from sam, this should probably just be passed as an option to bwa mem
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam

# sort and index for filtering
samtools sort -@ ${THREADS} ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
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