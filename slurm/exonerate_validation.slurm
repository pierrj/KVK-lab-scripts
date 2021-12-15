#!/bin/bash
#SBATCH --job-name=exonerate_validation
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/moryzae_pav/exonerate_validation

module load exonerate

LOST_GENOME=CH0043.fasta.masked
LOST_OG=OG0000000_protein.fasta

# exonerate --fsmmemory 2048 --softmasktarget --model protein2genome --showtargetgff yes --bestn 1 --showalignment no --verbose 0 --showvulgar no ${LOST_OG} ${LOST_GENOME} > missing.gff


tblastn -query ${LOST_OG} -subject ${LOST_GENOME} -out tblastn_test