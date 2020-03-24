#!/bin/bash
#SBATCH --job-name=samtools_sort
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/guy11_and_spikeins
samtools sort mergedandpe.CT_1A_bwamem.bam > sorted.mergedandpe.CT_1A_bwamem.bam
samtools index sorted.mergedandpe.CT_1A_bwamem.bam