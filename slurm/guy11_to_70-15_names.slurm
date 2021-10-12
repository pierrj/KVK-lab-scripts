#!/bin/bash
#SBATCH --job-name=guy11_to_70-15_names
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/gene_analysis/gene_names

query=guy11_fungap_out_12_28_20_prot.faa

blastp -query ${query} -db uniport_moryzae_70-15.fasta \
    -out ${query}.out -num_threads 20 -outfmt "6 qacc sacc evalue qstart qend sstart send qlen slen pident length stitle " -max_target_seqs 1


#80% sequence identity over 80% reciprocal sequence length

# for reciprocal, take length of alignment and make sure it is greater than 80% of both

awk '$10 >= 80 && $11/$8 >= 0.8 && $11/$9 >= 0.8' ${query}.out > ${query}.out.filtered

awk -v OFS='\t' '{print $1, substr($(NF-2),4)}' ${query}.out.filtered > ${query}.out.filtered.names

awk -v OFS='\t' '{print $1}' ${query}.out.filtered > annotated_genes

grep '>' $query | awk '{print substr($1, 2)}' > all_genes

cat annotated_genes all_genes | awk -v OFS='\t' '{ print $1, "unknown"}' | sort | uniq -u | cat - ${query}.out.filtered.names | sort -k1,1 | awk -v OFS='\t' '{print $2, $1}' > guy11_70_15_gene_names.txt