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
sample=IF_1A
samtools view -f 4 mergedandpe.${sample}_bwamem.bam | shuf -n 100 > test1.fasta
blastn -query test1.fasta -db /global/scratch/users/pierrj/blast/blastdb_ver4_curl/blast/db/v4/nt -out test1 -num_threads 20 -outfmt "6 qacc sacc evalue qstart qend sstart send stitle " -max_target_seqs 1
samtools view -f 4 sorted.mergedandpe.${sample}_bwamem.bam | shuf -n 100 > test2.fasta
blastn -query test2.fasta -db /global/scratch/users/pierrj/blast/blastdb_ver4_curl/blast/db/v4/nt -out test2 -num_threads 20 -outfmt "6 qacc sacc evalue qstart qend sstart send stitle " -max_target_seqs 1
samtools view -f 4 mergedandpe.${sample}_bwamem.bam | head -100 > test3.fasta
blastn -query test3.fasta -db /global/scratch/users/pierrj/blast/blastdb_ver4_curl/blast/db/v4/nt -out test3 -num_threads 20 -outfmt "6 qacc sacc evalue qstart qend sstart send stitle " -max_target_seqs 1
samtools view -f 4 mergedandpe.${sample}_bwamem.bam | head -100 > test4.fasta
blastn -query test4.fasta -db /global/scratch/users/pierrj/blast/blastdb_ver4_curl/blast/db/v4/nt -out test4 -num_threads 20 -outfmt "6 qacc sacc evalue qstart qend sstart send stitle " -max_target_seqs 1