#!/bin/bash
#SBATCH --job-name=test_70_15_fungap_to_uniprot
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/PAV_SV/PAV/test_assing_70_15_orthogroups

# from FunGAP
query=fungap_out_prot.faa

#from Kyungyong
# awk -F "," -v OFS="\n" '{ if ( $1 == "7") {print ">"$3, $5}}' Magnaporthe_Oryza_Structure_prediction_and_clustering_metadata.zenodo.csv > 70_15.ARTs.fasta
# awk -F "," -v OFS="\n" '{ if ( $1 == "8") {print ">"$3, $5}}' Magnaporthe_Oryza_Structure_prediction_and_clustering_metadata.zenodo.csv > 70_15.MAX.fasta
db=70_15.ARTs.fasta

makeblastdb -dbtype prot -in ${db}

blastp -query ${query} -db ${db} \
    -out ${query}.out -num_threads 20 -outfmt "6 qacc sacc evalue qstart qend sstart send qlen slen pident length stitle " -max_target_seqs 1

aw

#80% sequence identity over 80% reciprocal sequence length
# for reciprocal, take length of alignment and make sure it is greater than 80% of both
awk '$10 >= 80 && $11/$8 >= 0.5 && $11/$9 >= 0.5' ${query}.out > ${query}.out.filtered

awk -v OFS='\t' '{print $1, substr($(NF-2),4)}' ${query}.out.filtered > ${query}.out.filtered.names