#!/bin/bash
#SBATCH --job-name=map_spikeins_pUC18
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/pUC18_pluseverythingelse
sample=CT_1A
/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/${sample}_R1.fastq -r /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/${sample}_R2.fastq -1 tmp.seqprep.${sample}_R1.fastq.gz -2 tmp.seqprep.${sample}_R2.fastq.gz -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s tmp.merged.${sample}.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o tmp.trimmed.seqprep.${sample}_R1.fastq -p tmp.trimmed.seqprep.${sample}_R2.fastq tmp.seqprep.${sample}_R1.fastq.gz tmp.seqprep.${sample}_R2.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o tmp.trimmed.merged.${sample} tmp.merged.${sample}.fastq.gz
bwa mem -t 24 /global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_with_spikeins_pUC18_bwa tmp.trimmed.seqprep.${sample}_R1.fastq tmp.trimmed.seqprep.${sample}_R2.fastq -o tmp.seqprep.trimmed.${sample}_bwamem.sam
bwa mem -t 24 /global/scratch/users/pierrj/references/guy11_genome_baoetal2017_with_70-15_mito_with_spikeins_pUC18_bwa tmp.trimmed.merged.${sample} -o tmp.merged.trimmed.${sample}_bwamem.sam
samtools view -S -b tmp.seqprep.trimmed.${sample}_bwamem.sam > tmp.seqprep.trimmed.${sample}_bwamem.bam
samtools view -S -b tmp.merged.trimmed.${sample}_bwamem.sam > tmp.merged.trimmed.${sample}_bwamem.bam
bamtools merge -in tmp.seqprep.trimmed.${sample}_bwamem.bam -in tmp.merged.trimmed.${sample}_bwamem.bam -out mergedandpe.${sample}_bwamem.bam
samtools sort mergedandpe.${sample}_bwamem.bam > sorted.mergedandpe.${sample}_bwamem.bam
samtools index sorted.mergedandpe.${sample}_bwamem.bam

samtools view -f 81 -F 4 sorted.mergedandpe.${sample}_bwamem.bam pUC18 > tmp.reverseread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 145 -F 4 sorted.mergedandpe.${sample}_bwamem.bam pUC18 > tmp.reverseread2.mergedandpe.${sample}_bwamem.sam
samtools view -f 65 -F 20 sorted.mergedandpe.${sample}_bwamem.bam pUC18 > tmp.forwardread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 129 -F 20 sorted.mergedandpe.${sample}_bwamem.bam pUC18 > tmp.forwardread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.reverseread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.reverseread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.forwardread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread2.mergedandpe.${sample}_bwamem.sam tmp.forwardread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam

awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 19 ) || (b !~ /[DMIHS]/ && int(b) > 19)) print $0}' samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam > qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam > exactlytwice.qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam
samtools view -b -h <(cat <(samtools view -H mergedandpe.${sample}_bwamem.bam) exactlytwice.qualityfiltered.samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam) > SRs.mergedandpe.${sample}_bwamem.bam
bedtools bamtobed -i SRs.mergedandpe.${sample}_bwamem.bam | sort -k 4,4 -k 2,2 > SRs.mergedandpe.${sample}_bwamem.bed

rm tmp.*