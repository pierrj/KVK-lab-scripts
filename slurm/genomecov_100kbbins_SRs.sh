#!/bin/bash
#SBATCH --job-name=genomecov_100kbbins_SRs
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/codingregions_comparison/
sample=G3_1A
file="/global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/G3_1A/sorted.mergedandpe.${sample}_bwamem.bam"
samtools view -f 81 -F 4 sorted.mergedandpe.${sample}_bwamem.bam > tmp.reverseread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 145 -F 4 sorted.mergedandpe.${sample}_bwamem.bam > tmp.reverseread2.mergedandpe.${sample}_bwamem.sam
samtools view -f 65 -F 20 sorted.mergedandpe.${sample}_bwamem.bam > tmp.forwardread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 129 -F 20 sorted.mergedandpe.${sample}_bwamem.bam > tmp.forwardread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.reverseread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.reverseread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.forwardread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread2.mergedandpe.${sample}_bwamem.sam tmp.forwardread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam
