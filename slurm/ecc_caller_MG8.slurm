#!/bin/bash
#SBATCH --job-name=ecc_caller_MG&
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/CT_1A_MG8
sample=CT_1A
# /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/SeqPrep/SeqPrep/SeqPrep -f ${sample}_R1.fastq -r ${sample}_R2.fastq -1 tmp.seqprep.${sample}_R1.fastq.gz -2 tmp.seqprep.${sample}_R2.fastq.gz -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -s tmp.merged.${sample}.fastq.gz
# cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o tmp.trimmed.seqprep.${sample}_R1.fastq -p tmp.trimmed.seqprep.${sample}_R2.fastq tmp.seqprep.${sample}_R1.fastq.gz tmp.seqprep.${sample}_R2.fastq.gz
# cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o tmp.trimmed.merged.${sample} tmp.merged.${sample}.fastq.gz
# bwa mem -t 24 /global/scratch/users/pierrj/references/Magnaporthe_oryzae.MG8.dna.toplevel.noMT_bwa tmp.trimmed.seqprep.${sample}_R1.fastq tmp.trimmed.seqprep.${sample}_R2.fastq -o tmp.seqprep.trimmed.${sample}_bwamem.sam
# bwa mem -t 24 /global/scratch/users/pierrj/references/Magnaporthe_oryzae.MG8.dna.toplevel.noMT_bwa tmp.trimmed.merged.${sample} -o tmp.merged.trimmed.${sample}_bwamem.sam
# samtools view -S -b tmp.seqprep.trimmed.${sample}_bwamem.sam > tmp.seqprep.trimmed.${sample}_bwamem.bam
# samtools view -S -b tmp.merged.trimmed.${sample}_bwamem.sam > tmp.merged.trimmed.${sample}_bwamem.bam
# bamtools merge -in tmp.seqprep.trimmed.${sample}_bwamem.bam -in tmp.merged.trimmed.${sample}_bwamem.bam -out MG8.mergedandpe.${sample}_bwamem.bam
# samtools sort MG8.mergedandpe.${sample}_bwamem.bam > sorted.MG8.mergedandpe.${sample}_bwamem.bam
# samtools index sorted.MG8.mergedandpe.${sample}_bwamem.bam

# bedtools genomecov -d -ibam sorted.MG8.mergedandpe.${sample}_bwamem.bam > genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam

# cp /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/MG8.mergedandpe.${sample}_bwamem.bam /global/scratch/users/pierrj/eccDNA/pipeline_tests/CT_1A_MG8/mergedandpe.${sample}_bwamem.bam

# cp /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/sorted.MG8.mergedandpe.${sample}_bwamem.bam /global/scratch/users/pierrj/eccDNA/pipeline_tests/CT_1A_MG8/sorted.mergedandpe.${sample}_bwamem.bam

# cp /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/sorted.MG8.mergedandpe.${sample}_bwamem.bam.bai /global/scratch/users/pierrj/eccDNA/pipeline_tests/CT_1A_MG8/sorted.mergedandpe.${sample}_bwamem.bam.bai

# cp /global/scratch/users/pierrj/eccDNA/stress_experiments/mag_stress/CT_1A/genomecoverage.MG8.mergedandpe.${sample}_bwamem.bam /global/scratch/users/pierrj/eccDNA/pipeline_tests/CT_1A_MG8/genomecoverage.mergedandpe.${sample}_bwamem.bam

samtools view -b -q 1 sorted.mergedandpe.${sample}_bwamem.bam > guy11.sorted.mergedandpe.${sample}_bwamem.bam
samtools view -f 81 -F 4 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.reverseread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 145 -F 4 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.reverseread2.mergedandpe.${sample}_bwamem.sam
samtools view -f 65 -F 20 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.forwardread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 129 -F 20 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.forwardread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.reverseread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.reverseread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.forwardread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread2.mergedandpe.${sample}_bwamem.sam tmp.forwardread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam

samtools view guy11.sorted.mergedandpe.${sample}_bwamem.bam | awk '{ if (($2 == 83 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 99 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H mergedandpe.${sample}_bwamem.bam) - | samtools view -b -h - > tmp.outwardfacing.mergedandpe.${sample}_bwamem.bam
bedtools bamtobed -i tmp.outwardfacing.mergedandpe.${sample}_bwamem.bam | sort -k 4,4 > tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed
mv tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old
awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,substr($4, 1, length($4)-2),$5,$6}' tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old > tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old.trimmed
awk 'NR==FNR{a[$4]++; next} a[$4]==2' tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old.trimmed tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old.trimmed > outwardfacing.mergedandpe.${sample}_bwamem.bed

ipcluster start -n $SLURM_NTASKS --cluster-id="slurm-${SLURM_JOBID}" &
sleep 300
ipython /global/home/users/pierrj/git/python/ecc_callerv3_flexiblescaffolds.py samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam outwardfacing.mergedandpe.${sample}_bwamem.bed genomecoverage.mergedandpe.${sample}_bwamem.bed ${sample} 7 '' 0 1 ''
ipcluster stop --cluster-id="slurm-${SLURM_JOBID}"