#!/bin/bash
#SBATCH --job-name=WGS-allSRs
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
cd /global/scratch/users/pierrj/genome/2018_foster/
sample=ERR2660591
igvtools count mergedandpe.${sample}_bwamem.bam /global/scratch/users/pierrj/eccDNA/pipeline_tests/peaks_gone/WGS.illumina.bam.tdf /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta
# samtools view -b -q 1 mergedandpe.${sample}_bwamem.bam > unique.mergedandpe.${sample}_bwamem.bam
igvtools count unique.mergedandpe.${sample}_bwamem.bam /global/scratch/users/pierrj/eccDNA/pipeline_tests/peaks_gone/WGS.illumina.unique.bam.tdf /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta

# samtools view -f 81 -F 4 mergedandpe.${sample}_bwamem.bam > tmp.reverseread1.mergedandpe.${sample}_bwamem.sam
# samtools view -f 145 -F 4 mergedandpe.${sample}_bwamem.bam > tmp.reverseread2.mergedandpe.${sample}_bwamem.sam
# samtools view -f 65 -F 20 mergedandpe.${sample}_bwamem.bam > tmp.forwardread1.mergedandpe.${sample}_bwamem.sam
# samtools view -f 129 -F 20 mergedandpe.${sample}_bwamem.bam > tmp.forwardread2.mergedandpe.${sample}_bwamem.sam
# awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.reverseread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
# awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.reverseread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
# awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.forwardread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
# awk 'NR==FNR{a[$1]++; next} a[$1]==2' tmp.forwardread2.mergedandpe.${sample}_bwamem.sam tmp.forwardread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
# cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > exactlytwice.all.mergedandpe.${sample}_bwamem.bam

# samtools view -b <(cat <(samtools view -H mergedandpe.${sample}_bwamem.bam) exactlytwice.all.mergedandpe.${sample}_bwamem.bam ) > exactlytwice.all.mergedandpe.${sample}_bwamem.actual.bam
# bedtools bamtobed -i exactlytwice.all.mergedandpe.${sample}_bwamem.actual.bam | sort -k 1,1 -k2n > exactlytwice.all.mergedandpe.${sample}_bwamem.bed

# igvtools count exactlytwice.all.mergedandpe.${sample}_bwamem.bed /global/scratch/users/pierrj/eccDNA/pipeline_tests/peaks_gone/WGS.illumina.allSRs.bed.tdf /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta


# samtools view -b <(cat <(samtools view -H mergedandpe.${sample}_bwamem.bam) samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam ) > samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.actual.bam
# bedtools bamtobed -i samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.actual.bam | sort -k 1,1 -k2n > samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.actual.bed

# igvtools count samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.actual.bed /global/scratch/users/pierrj/eccDNA/pipeline_tests/peaks_gone/WGS.illumina.samechromSRs.bed.tdf /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta