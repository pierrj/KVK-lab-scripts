#!/bin/bash
#SBATCH --job-name=validate_missing_genes
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
cd /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/

while read genome; do
blastn -query ${genome}_missing_genes.fasta -db genomes/${genome}.blastdb \
    -num_threads ${SLURM_NTASKS} -outfmt "6 qacc sacc evalue qlen qstart qend sstart send  stitle " \
    -max_target_seqs 1 | awk '{if (($6-$5)/$4>0.25) {print $1}}' | sort | uniq > ${genome}.not_missing_genes
done < moryzae_oryza_genomes_mapfile