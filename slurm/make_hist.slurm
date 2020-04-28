#!/bin/bash
#SBATCH --job-name=make_hist
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/distribution
awk '{print $3-$2}' ../merging/parallel.confirmed.G3_1.bed > distribution.G3_1.unsorted
R --no-save
chol <- read.table("distribution.G3_1.unsorted", header = FALSE)
tiff("distribution.G3_1.unsorted.tiff", units="in", width=5, height=5, res=300)
hist(chol$V1, breaks=1000)
dev.off()