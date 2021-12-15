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

LOST_GENOME=CH0043.fasta
LOST_OG=OG0000323_protein.fasta
LOST_OG=OG0000000_protein.fasta

exonerate --model protein2genome --showtargetgff yes --bestn 1 --showalignment no --verbose 0 --showvulgar no ${LOST_OG} ${LOST_GENOME} > missing.gff


# tblastn -query ${LOST_OG} -subject ${LOST_GENOME} -out tblastn_test_positive_control

tblastn -query ${LOST_OG} -subject ${LOST_GENOME} \
    -outfmt "6 qacc sacc evalue qlen qstart qend sstart send  stitle " \
    -max_target_seqs 1 > tblastn_test_positive_control_tabs


tblastn -query ${LOST_OG} -subject ${LOST_GENOME} \
    -evalue 1e-10 -max_intron_length 3000 \
    -outfmt "6 qacc sacc evalue qlen qstart qend sstart send  pident" \
    -max_target_seqs 1 > tblastn_test_positive

tblastn -query ${LOST_OG} -subject ${LOST_GENOME} \
    -evalue 1e-10 -max_intron_length 3000 \
    -outfmt "6 qacc sacc evalue qlen qstart qend sstart send  pident qcovs qcovhsp"  \
    -max_target_seqs 1 > tblastn_test_positive