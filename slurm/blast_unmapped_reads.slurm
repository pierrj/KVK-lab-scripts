#!/bin/bash
#SBATCH --job-name=blast_unmapped_reads
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
cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/infected_and_spikeins
blastn -query  unmapped100.sorted.mergedandpe.IF_1A_bwamem.fasta -db nt -out unmapped100.IF_1A.blastn -num_threads 20 -outfmt 6 -max_target_seqs 1