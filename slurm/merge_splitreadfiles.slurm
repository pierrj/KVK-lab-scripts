#!/bin/bash
#SBATCH --job-name=ecc_caller_test
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
cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/rawdata/illumina/G3_1A_tests
samtools view -f 81 -F 4 mergedandpe.G3_1A_bwamem.bam > reverseread1.mergedandpe.G3_1A_bwamem.bam
samtools view -f 145 -F 4 mergedandpe.G3_1A_bwamem.bam > reverseread2.mergedandpe.G3_1A_bwamem.bam
samtools view -f 65 -F 20 mergedandpe.G3_1A_bwamem.bam > forwardread1.mergedandpe.G3_1A_bwamem.bam
samtools view -f 129 -F 20 mergedandpe.G3_1A_bwamem.bam > forwardread2.mergedandpe.G3_1A_bwamem.bam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' reverseread1.mergedandpe.G3_1A_bwamem.bam reverseread1.mergedandpe.G3_1A_bwamem.bam > samechromosome.exactlytwice.reverseread1.mergedandpe.G3_1A_bwamem.bam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' reverseread2.mergedandpe.G3_1A_bwamem.bam reverseread2.mergedandpe.G3_1A_bwamem.bam > samechromosome.exactlytwice.reverseread2.mergedandpe.G3_1A_bwamem.bam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' forwardread1.mergedandpe.G3_1A_bwamem.bam forwardread1.mergedandpe.G3_1A_bwamem.bam > samechromosome.exactlytwice.forwardread1.mergedandpe.G3_1A_bwamem.bam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' forwardread2.mergedandpe.G3_1A_bwamem.bam forwardread2.mergedandpe.G3_1A_bwamem.bam > samechromosome.exactlytwice.forwardread2.mergedandpe.G3_1A_bwamem.bam
cat samechromosome.exactlytwice.reverseread1.mergedandpe.G3_1A_bwamem.bam samechromosome.exactlytwice.reverseread2.mergedandpe.G3_1A_bwamem.bam samechromosome.exactlytwice.forwardread1.mergedandpe.G3_1A_bwamem.bam samechromosome.exactlytwice.forwardread2.mergedandpe.G3_1A_bwamem.bam > samechromosome.exactlytwice.all.mergedandpe.G3_1A_bwamem.bam