#!/bin/bash
#SBATCH --job-name=singleend_poption_test
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
SAMPLE="CT_1A"
READONE="/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/CT_1A_R1.fastq"
READTWO="/global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/CT_1A_R2.fastq"
GENOME_DB="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_with_spikeins_tweaked_bwa"
THREADS="24"
MAPFILE="/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames"

# cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/single_end_test
# gzip -c ${READONE} > tmp.seqprep.${SAMPLE}_R1.fastq.gz
# gzip -c ${READTWO} > tmp.seqprep.${SAMPLE}_R2.fastq.gz
# cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq tmp.seqprep.${SAMPLE}_R1.fastq.gz tmp.seqprep.${SAMPLE}_R2.fastq.gz
# cat tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq > tmp.trimmed.seqprep.${SAMPLE}.fastq
# bwa mem -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
# samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam
# samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
# samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

# samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam

# rm tmp*

# cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/p_option_map_test
# gzip -c ${READONE} > tmp.seqprep.${SAMPLE}_R1.fastq.gz
# gzip -c ${READTWO} > tmp.seqprep.${SAMPLE}_R2.fastq.gz
# cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq tmp.seqprep.${SAMPLE}_R1.fastq.gz tmp.seqprep.${SAMPLE}_R2.fastq.gz
# bwa mem -P -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
# samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam
# samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
# samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

# samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam

# rm tmp*

# mkdir /global/scratch/users/pierrj/eccDNA/pipeline_tests/a_option_map_test
# cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/a_option_map_test
# gzip -c ${READONE} > tmp.seqprep.${SAMPLE}_R1.fastq.gz
# gzip -c ${READTWO} > tmp.seqprep.${SAMPLE}_R2.fastq.gz
# cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq tmp.seqprep.${SAMPLE}_R1.fastq.gz tmp.seqprep.${SAMPLE}_R2.fastq.gz
# cat tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq > tmp.trimmed.seqprep.${SAMPLE}.fastq
# bwa mem -a -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
# samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam
# samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
# samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

# samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam

# rm tmp*

mkdir /global/scratch/users/pierrj/eccDNA/pipeline_tests/star_map_test
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/star_map_test
gzip -c ${READONE} > tmp.seqprep.${SAMPLE}_R1.fastq.gz
gzip -c ${READTWO} > tmp.seqprep.${SAMPLE}_R2.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq tmp.seqprep.${SAMPLE}_R1.fastq.gz tmp.seqprep.${SAMPLE}_R2.fastq.gz
STAR --runThreadN 24 \
    --genomeDir /global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_with_spikeins_tweaked_star  \
    --readFilesIn tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq \
    --outSAMtype BAM SortedByCoordinate --runMode alignReads --outFilterMultimapNmax 1000 \
    --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random \
    --winAnchorMultimapNmax 1000 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 350

rm tmp*

mdkir /global/scratch/users/pierrj/eccDNA/pipeline_tests/bwamem_newoptions
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/bwamem_newoptions
gzip -c ${READONE} > tmp.seqprep.${SAMPLE}_R1.fastq.gz
gzip -c ${READTWO} > tmp.seqprep.${SAMPLE}_R2.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq tmp.seqprep.${SAMPLE}_R1.fastq.gz tmp.seqprep.${SAMPLE}_R2.fastq.gz
cat tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq > tmp.trimmed.seqprep.${SAMPLE}.fastq
bwa mem -t ${THREADS} -c 50000 -T 20 ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam

rm tmp*