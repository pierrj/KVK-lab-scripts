#!/bin/bash
#SBATCH --job-name=callingwithmulti-mappingreads-allSRs
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/peaks_gone/ecc_caller_with_repeatmapping
sample=G3_1A
samtools view -f 81 -F 4 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.reverseread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 145 -F 4 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.reverseread2.mergedandpe.${sample}_bwamem.sam
samtools view -f 65 -F 20 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.forwardread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 129 -F 20 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.forwardread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.reverseread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.reverseread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.forwardread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.forwardread2.mergedandpe.${sample}_bwamem.sam tmp.forwardread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > exactlytwice.all.mergedandpe.${sample}_bwamem.bam

samtools view -b <(cat <(samtools view -H mergedandpe.${sample}_bwamem.bam) exactlytwice.all.mergedandpe.${sample}_bwamem.bam ) > exactlytwice.all.mergedandpe.${sample}_bwamem.actual.bam
bedtools bamtobed -i exactlytwice.all.mergedandpe.${sample}_bwamem.actual.bam | sort -k 4,4 -k 2,2 > exactlytwice.all.mergedandpe.${sample}_bwamem.bed

igvtools count exactlytwice.all.mergedandpe.${sample}_bwamem.bed /global/scratch/users/pierrj/eccDNA/pipeline_tests/peaks_gone/illumina.allSRs.G3_1A.bed.tdf /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta