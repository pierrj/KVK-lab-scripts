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
cd /global/scratch/users/pierrj/eccDNA/2018_moller/SRR6315396
sample=SRR6315396
cutadapt --interleaved -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o trimmed.${sample}_R1.fastq -p trimmed.${sample}_R2.fastq ${sample}_1.fastq ${sample}_2.fastq
/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f trimmed.${sample}_R1.fastq -r trimmed.${sample}_R2.fastq -1 seqprep.trimmed.${sample}_R1.fastq -2 seqprep.trimmed.${sample}_R2.fastq -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s merged.trimmed.${sample}
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.bwa seqprep.trimmed.${sample}_R1.fastq seqprep.trimmed.${sample}_R2.fastq -o seqprep.trimmed.${sample}_bwamem.sam
bwa mem -t 24 /global/scratch/users/pierrj/eccDNA/2018_moller/references/GCF_000001405.25_GRCh37.p13_genomic.bwa merged.trimmed.${sample} -o merged.trimmed.${sample}_bwamem.sam
samtools view -S -b seqprep.trimmed.${sample}_bwamem.sam > seqprep.trimmed.${sample}_bwamem.bam
samtools view -S -b merged.trimmed.${sample}_bwamem.sam > merged.trimmed.${sample}_bwamem.bam
bamtools merge -in seqprep.trimmed.${sample}_bwamem.bam -in merged.trimmed.${sample}_bwamem.bam -out mergedandpe.${sample}_bwamem.bam

## get splitread input file
samtools view -f 81 -F 4 mergedandpe.${sample}_bwamem.bam > reverseread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 145 -F 4 mergedandpe.${sample}_bwamem.bam > reverseread2.mergedandpe.${sample}_bwamem.sam
samtools view -f 65 -F 20 mergedandpe.${sample}_bwamem.bam > forwardread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 129 -F 20 mergedandpe.${sample}_bwamem.bam > forwardread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' reverseread1.mergedandpe.${sample}_bwamem.sam reverseread1.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' reverseread2.mergedandpe.${sample}_bwamem.sam reverseread2.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' forwardread1.mergedandpe.${sample}_bwamem.sam forwardread1.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' forwardread2.mergedandpe.${sample}_bwamem.sam forwardread2.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
cat  samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam

samtools view mergedandpe.${sample}_bwamem.bam | awk '{ if (($2 == 83 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 99 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H mergedandpe.${sample}_bwamem.bam) - | samtools view -b -h - > outwardfacing.mergedandpe.${sample}_bwamem.bam
bedtools bamtobed -i outwardfacing.mergedandpe.${sample}_bwamem.bam | sort -k 4,4 > outwardfacing.mergedandpe.${sample}_bwamem.bed

samtools sort mergedandpe.${sample}_bwamem.bam -o sorted.mergedandpe.${sample}_bwamem.bam
bedtools genomecov -d -ibam sorted.mergedandpe.${sample}_bwamem.bam > genomecoverage.mergedandpe.${sample}_bwamem.bed

ipcluster start -n $SLURM_NTASKS --cluster-id="slurm-${SLURM_JOBID}" &
sleep 45
ipython /global/home/users/pierrj/git/python/ecc_caller.py samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam outwardfacing.mergedandpe.${sample}_bwamem.bed genomecoverage.mergedandpe.${sample}_bwamem.bed ${sample}
ipcluster stop --cluster-id="slurm-${SLURM_JOBID}"